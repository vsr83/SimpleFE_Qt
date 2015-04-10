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

#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QStatusBar>
#include <QMainWindow>
#include <QMenu>
#include <QMenuBar>
#include <QLabel>
#include <QAction>
#include <QSet>
#include <QSplitter>
#include <map>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QSpinBox>

#include "physlist.h"
#include "meshplot.h"
#include "mesh.h"
#include "region.h"
#include "partition.h"

class MainWindow : public QMainWindow {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = 0);
    ~MainWindow();

public slots:
    void updateSolution();
    void updateContours(int _numContours);
protected:
    void keyPressEvent(QKeyEvent *);
private slots:
    void about();
    void setAntiAlias(bool aa);
    void setAxisEqual(bool ae);
    void setColorbar (bool sc);
    void setShowMesh (bool sm);
    void updateSol();
    void openMesh();
    void openSolution();
    void hilightDialog();
    void physChanged(std::map<int, Region> _physmap, int phys);
    void physSelected(int phys);

    void undo(std::map <int, Region> _physmap);
    void redo(std::map <int, Region> _physmap);
    void reset(std::map <int, Region> _physmap);
private:
    QLabel *statusLabel, *statusLabel2;
    QMenu *fileMenu;
    QMenu *toolsMenu;
    QMenu *optionsMenu;
    QMenu *helpMenu;
    QAction *openMeshAction, *openSolutionAction, *updateSolutionAction;
    QAction *quitAction;
    QAction *antiaAction, *axiseqAction, *colorbarAction, *hilightAction, *showmeshAction;
    QAction *aboutAction;
    QSplitter *splitter, *hsplitter;
    QWidget *hboxWidget, *vboxWidget;
    QHBoxLayout *hboxLayout;
    QVBoxLayout *vboxLayout;
    MeshPlot *meshplot;
    physList *physlist;
    QLabel *contourLabel;
    QSpinBox *contourSpinBox;

    MeshFile *meshfile;
    Mesh     *mesh;
    Partition *meshpart;

    int numContours;

    QMap <int, QColor> hilightMap;
    QVector <double> solution;
    int numDoF;
    std::map <int, int>    partition;
    std::map <int, Region> physmap;
    double mu0;
};

#endif // MAINWINDOW_H
