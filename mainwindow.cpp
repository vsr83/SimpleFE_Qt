/* SimpleFE_Qt - A Simple FE Solver with a graphical interface.
   Copyright (C) 2015 Ville Räisänen <vsr at vsr.name>

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "mainwindow.h"
#include <QDebug>
#include <QKeyEvent>
#include <QMessageBox>
#include <QFileDialog>
#include <QFileInfo>
#include <QInputDialog>
#include <math.h>

#include "assembly.h"
#include "partition.h"

MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent) {
    mu0         = 4*M_PI*1e-7;
    meshfile    = 0;
    mesh        = 0;
    meshpart    = 0;
    numDoF      = 0;
    numContours = 20;

    statusLabel = new QLabel("Status");
    statusLabel2 = new QLabel("--");

    statusBar()->addWidget(statusLabel, 1);
    statusBar()->addWidget(statusLabel2);

    openMeshAction = new QAction(tr("Open &Mesh"), this);
    openMeshAction->setShortcut(tr("Ctrl+M"));
    connect(openMeshAction, SIGNAL(triggered()), SLOT(openMesh()));

    openSolutionAction = new QAction(tr("Open &Solution"), this);
    openSolutionAction->setShortcut(tr("Ctrl+S"));
    connect(openSolutionAction, SIGNAL(triggered()), SLOT(openSolution()));

    updateSolutionAction = new QAction(tr("&Update Solution"), this);
    updateSolutionAction->setShortcut(tr("Ctrl+U"));
    connect(updateSolutionAction, SIGNAL(triggered()), SLOT(updateSol()));


    quitAction = new QAction(tr("&Quit"), this);
    quitAction->setShortcut(QKeySequence::Quit);
    quitAction->setShortcut(tr("Ctrl+Q"));
    connect(quitAction, SIGNAL(triggered()), this, SLOT(close()));

    fileMenu = menuBar()->addMenu(tr("&File"));
    fileMenu->addAction(openMeshAction);
    fileMenu->addAction(openSolutionAction);
    fileMenu->addAction(updateSolutionAction);
    fileMenu->addAction(quitAction);
//    toolsMenu = menuBar()->addMenu(tr("&Tools"));
    optionsMenu = menuBar()->addMenu(tr("&Options"));

    showmeshAction = new QAction(tr("Show Mesh"), this);
    showmeshAction->setCheckable(true);
    showmeshAction->setChecked(true);
    optionsMenu->addAction(showmeshAction);
    connect(showmeshAction, SIGNAL(toggled(bool)), this, SLOT(setShowMesh(bool)));

    antiaAction = new QAction(tr("&Antialiasing"), this);
    antiaAction->setCheckable(true);
    antiaAction->setChecked(true);
    optionsMenu->addAction(antiaAction);
    connect(antiaAction, SIGNAL(toggled(bool)), this, SLOT(setAntiAlias(bool)));

    axiseqAction = new QAction(tr("Force Axis Equal"), this);
    axiseqAction->setCheckable(true);
    axiseqAction->setChecked(true);
    optionsMenu->addAction(axiseqAction);
    connect(axiseqAction, SIGNAL(toggled(bool)), this, SLOT(setAxisEqual(bool)));

    colorbarAction = new QAction(tr("Show Colorbar"), this);
    colorbarAction->setCheckable(true);
    colorbarAction->setChecked(true);
    optionsMenu->addAction(colorbarAction);
    connect(colorbarAction, SIGNAL(toggled(bool)), this, SLOT(setColorbar(bool)));

    hilightAction = new QAction(tr("Hilight Region"), this);
    hilightAction->setShortcut(tr("Ctrl+H"));
    optionsMenu->addAction(hilightAction);
    connect(hilightAction, SIGNAL(triggered()), SLOT(hilightDialog()));

    menuBar()->addSeparator();
    helpMenu = menuBar()->addMenu(tr("&Help"));
    aboutAction = new QAction(tr("About"), this);
    helpMenu->addAction(aboutAction);
    connect(aboutAction, SIGNAL(triggered()), this, SLOT(about()));

    meshplot = new MeshPlot;
    meshplot->setFocus();
    meshplot->setFocusPolicy(Qt::StrongFocus);

    std::map <int, Region> physmap;
    Region empty;
    empty.name = std::string("Empty");
    empty.nu = 0;
    empty.Js = 0;

    physmap[1] = empty;

    physlist = new physList(physmap);
    connect(physlist, SIGNAL(physChanged(std::map<int,Region>,int)),
            this,       SLOT(physChanged(std::map<int,Region>,int)));
    connect(physlist, SIGNAL(physSelected(int)),
            this,       SLOT(physSelected(int)));
    connect(physlist, SIGNAL(undoClicked(std::map<int,Region>)),
            this,       SLOT(undo(std::map <int, Region>)));
    connect(physlist, SIGNAL(redoClicked(std::map <int, Region>)),
            this,       SLOT(redo(std::map <int, Region>)));
    connect(physlist, SIGNAL(resetClicked(std::map <int, Region>)),
            this,       SLOT(reset(std::map <int, Region>)));

    splitter = new QSplitter(Qt::Vertical);
    hboxLayout = new QHBoxLayout;
    vboxLayout = new QVBoxLayout;
    hboxWidget = new QWidget;
    vboxWidget = new QWidget;
    hboxWidget->setLayout(hboxLayout);
    vboxWidget->setLayout(vboxLayout);

    hboxLayout->addWidget(physlist);
    hboxLayout->addWidget(vboxWidget);
    hboxLayout->setContentsMargins(0, 0, 0, 0);
    splitter->addWidget(meshplot);
    splitter->addWidget(hboxWidget);

    contourLabel = new QLabel(tr("numContours"));
    vboxLayout->addWidget(contourLabel);
    contourSpinBox = new QSpinBox;
    contourSpinBox->setRange(1, 1000);
    contourSpinBox->setValue(numContours);
    vboxLayout->addWidget(contourSpinBox);
    vboxLayout->addStretch();

    connect(contourSpinBox, SIGNAL(valueChanged(int)),
            this,             SLOT(updateContours(int)));

    physlist->show();
    this->setCentralWidget(splitter);
    resize(800, 600);
}

void
MainWindow::undo(std::map <int, Region> _physmap) {
    physmap = _physmap;
    updateSol();
}

void
MainWindow::redo(std::map <int, Region> _physmap) {
    physmap = _physmap;
    updateSol();
}

void
MainWindow::reset(std::map <int, Region> _physmap) {
    physmap = _physmap;
    updateSol();
}


void
MainWindow::physSelected(int phys) {
    qDebug() << "MAIN - physSelected";
}

void
MainWindow::updateSol() {
    updateSolution();
    updateContours(numContours);
}
void
MainWindow::updateSolution() {
    qDebug() << "updateSolution()";
    /*
    for (std::map<int, Region>::iterator it = physmap.begin();
         it != physmap.end(); ++it) {
        int phys = it->first;
        Region reg = it->second;

        qDebug() << phys << reg.name.data() << reg.nu << reg.Js;
    }
    */

    // Assemble the global stiffness and mass matrices.
    Assembly ass(mesh, 6, physmap);
    int num_free = meshpart->parts_size[0];
    int num_bnd  = meshpart->parts_size[1];

    Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
    Eigen::SparseMatrix<double> sol_sparse, sol_sparsefull;
    Eigen::MatrixXd             sol_dense,  sol_densefull;
    Eigen::SparseMatrix<double> S_FI, S_FF, F_F;

    // Partition nodes into free and bounded nodes (see partition.h).
    S_FI = Eigen::SparseMatrix<double>(num_free, num_bnd);
    S_FF = Eigen::SparseMatrix<double>(num_free, num_free);
    F_F  = Eigen::SparseMatrix<double>(num_free, 1);
    S_FF = meshpart->part_left(0, 0)*(*ass.global_stiff)*meshpart->part_right(0, 0);
    S_FI = meshpart->part_left(0, 1)*(*ass.global_stiff)*meshpart->part_right(0, 1);
    F_F  = meshpart->part_left(0, 0) * ass.excitation;

    // Solve the system of equations.
    solver.compute(S_FF);
    if (solver.info() != Eigen::Success) {
        QMessageBox::critical(this, tr("Eigen::SparseLU!"),
                              tr("SpraseLU Decomposition Failed!\n"
                                 "."), QMessageBox::Ok, 0);
        std::cerr << "SparseLU decomposition failed." << std::endl;
        return;
    }

    sol_sparse    = solver.solve(F_F);
    sol_sparsefull= (meshpart->part_right(0, 0)) * sol_sparse;
    sol_dense     = Eigen::MatrixXd(sol_sparse);
    sol_densefull = Eigen::MatrixXd(sol_sparsefull);

    //std::cout << sol_sparse;
    if (solver.info() != Eigen::Success) {
        QMessageBox::critical(this, tr("Eigen::SparseLU!"),
                              tr("SpraseLU Solution Failed!\n"
                                 "."), QMessageBox::Ok, 0);
        std::cerr << "Solution failed." << std::endl;
        return;
    }

    numDoF = mesh->num_nodes;
    solution.clear();
    solution.reserve(numDoF);
    for (int ind_dof = 0; ind_dof < numDoF; ind_dof++) {
       solution.push_back(sol_densefull(ind_dof, 0));
    }
}

void
MainWindow::updateContours(int _numContours) {
    numContours = _numContours;

    // Remove contour layer(s) (layers with > 0 index):
    if (meshplot->visdata.numLayers >= 2) {
        for (int ind_layer = meshplot->visdata.numLayers-1; ind_layer > 0; ind_layer--) {
            delete meshplot->visdata.layers[ind_layer];
            meshplot->visdata.layers.remove(ind_layer);
        }
        meshplot->visdata.numLayers = 1;
    }

    Layer *contourLayer = new Layer;
    contourLayer->createContours(mesh, solution, numContours);
    meshplot->visdata.layers.push_back(contourLayer);
    meshplot->visdata.numLayers++;
    meshplot->refreshPixmap();
}

void
MainWindow::physChanged(std::map<int, Region> _physmap, int phys) {
    qDebug() << "MAIN - physChanged";

    physmap = _physmap;
    updateSolution();
    updateContours(numContours);
}

void
MainWindow::about() {
    QMessageBox::about(this, tr("About"),
                       tr("<h2>SimpleFE Postprocessor 0.1</h2>"
                          "<p>Copyright &copy; 2015 Ville Räisänen vsr@vsr.name</p>"
                          "<p>See COPYING for license information</p>"));
}

void
MainWindow::setShowMesh(bool sm) {
    if (meshplot->visdata.numLayers > 0) {
        meshplot->visdata.layers[0]->linesVisible = sm;
        meshplot->refreshPixmap();
    }
}

void
MainWindow::setAntiAlias(bool aa) {
    meshplot->setAntiAliasing(aa);
    meshplot->refreshPixmap();
}

void
MainWindow::setAxisEqual(bool ae) {
    meshplot->setAxisEqual(ae);
    if (ae) {
        meshplot->fixAspectRatio();
        meshplot->refreshPixmap();
    }
}

void
MainWindow::setColorbar(bool sc) {
    meshplot->setShowColorbar(sc);
    meshplot->refreshPixmap();
}

void
MainWindow::hilightDialog() {
    bool ok;
    int phys = QInputDialog::getInt(this, tr("Hilight Region"),
                                    tr("Hilight Region"),
                                    1, 1, 1e9, 1, &ok);
    if (ok && meshplot->visdata.numLayers > 0 && mesh) {
        QColor hilightColor(160, 160, 160);

        if (hilightMap.contains(phys)) {
            hilightMap.remove(phys);
        } else {
            hilightMap[phys] = hilightColor;
        }
        meshplot->visdata.layers[0]->polygons.clear();
        meshplot->visdata.layers[0]->numPolygons = 0;
        meshplot->visdata.layers[0]->colorPhysicals(mesh, hilightMap);
        meshplot->refreshPixmap();
    }
}

void
MainWindow::openMesh() {
    // Generate a dialog for the selection of the GMsh mesh file.
    QString filename = QFileDialog::getOpenFileName(this, tr("Open Mesh File"),
                                                    ".", tr("Mesh Files (*.msh)"));
    QFileInfo info(filename);
    if (!info.exists()) {
        return;
    }

    // Load the mesh file:
    if (mesh) {
        hilightMap.clear();
        delete meshfile;
        delete mesh;
        meshplot->clearData();
        meshfile = new MeshFile(qPrintable(filename));
        mesh = new Mesh(meshfile);
    } else {
        meshfile = new MeshFile(qPrintable(filename));
        mesh = new Mesh(meshfile);
    }

    // Update status bar.
    QString qs;
    QString file = QFileInfo(filename).fileName();
    qs = qs.sprintf("%s, %d Nodes, %d Triangles", qPrintable(file),
                          mesh->num_nodes, mesh->num_triangles);
    statusLabel->setText(qs);

    // Initialize the meshplot for the visualization of the mesh.
    meshplot->loadMesh(mesh);
    meshplot->show();
    meshplot->fixAspectRatio();
    meshplot->refreshPixmap();

    // Generate initial material parameters for the problem and partition the
    // physical regions to free and bounded nodes. In this code zero Dirichlet
    // condition is assumed on all physical lines.

    physmap.clear();
    for (unsigned int ind_tri=0; ind_tri < mesh->num_triangles; ind_tri++) {
        Mesh_Element *elem = mesh->elements[mesh->triangles[ind_tri]];
        Region reg(std::string("Empty"), 1/mu0, 0, 0);
        physmap[elem->physical] = reg;
        partition[elem->physical] = 0;
    }
    for (unsigned int ind_line=0; ind_line < mesh->num_lines; ind_line++) {
        Mesh_Element *elem = mesh->elements[mesh->lines[ind_line]];
        for (unsigned int ind_node=0; ind_node < elem->nnodes; ind_node++) {
            partition[elem->physical] = 1;
        }
    }

    // Generate the sparse matrices used for the partitioning of the global
    // stiffness matrix and the excitation vector. Change in the material
    // parameters or excitation does not influence these matrices.

    if (meshpart) delete meshpart;
    meshpart = new Partition(mesh, partition);

    // Update the list of physical data with initial data.
    physlist->updateList(physmap, 0);
    physlist->setFirst(physmap);
}

void
MainWindow::openSolution() {
    if (!mesh) {
        QMessageBox::critical(this, tr("No Mesh File Open!"),
                              tr("Cannot load solution\n"
                                 "No mesh file open."), QMessageBox::Ok, 0);
    } else {
        QString filename = QFileDialog::getOpenFileName(this, tr("Open Solution File"),
                                                    ".", tr("Solution Files *"));
        QFileInfo info(filename);
        if (!info.exists()) {
            return;
        }

        QFile file(filename);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
            return;
        double numDoF_tmp = 0;
        QVector <double> solution_tmp;
        solution_tmp.reserve(mesh->num_nodes);
        while (!file.atEnd()) {
            QByteArray line = file.readLine();
            line.remove(line.length()-1, 1);
            solution_tmp.push_back(line.toDouble());
            numDoF_tmp++;
        }
        bool ok;
        int numContours = QInputDialog::getInt(this, tr("Number of Contour Lines"),
                                               tr("Number of Contour Lines"),
                                               20, 1, 1000, 1, &ok);
        if (ok) {
            if (numDoF_tmp == mesh->num_nodes) {
                if (meshplot->visdata.numLayers >= 2) {
                    for (int ind_layer = meshplot->visdata.numLayers-1; ind_layer > 0; ind_layer++) {
                        delete meshplot->visdata.layers[ind_layer];
                        meshplot->visdata.layers.remove(ind_layer);
                    }
                    meshplot->visdata.numLayers = 1;
                }
                numDoF = numDoF_tmp;
                solution = solution_tmp;

                Layer *contourLayer = new Layer;
                contourLayer->createContours(mesh, solution, numContours);
                meshplot->visdata.layers.push_back(contourLayer);
                meshplot->visdata.numLayers++;
                meshplot->refreshPixmap();
            } else {
                QMessageBox::critical(this, tr("Invalid File!"),
                                      tr("Number of lines doesn't match with the\n"
                                         "number of mesh nodes."), QMessageBox::Ok, 0);
            }
        }
    }
}

void
MainWindow::keyPressEvent(QKeyEvent *event) {
    qDebug() << event->key();
}

MainWindow::~MainWindow() {
    if (meshfile) {
        delete meshfile;
        meshfile = 0;
    }
    if (mesh) {
        delete mesh;
        mesh = 0;
    }
}
