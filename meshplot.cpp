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


#include "meshplot.h"
#include "mesh.h"

#include <QDebug>
#include <QStylePainter>
#include <QWheelEvent>
#include <QMouseEvent>
#include <QPoint>
#include <QPointF>
#include <cmath>

MeshPlot::MeshPlot(QWidget *parent) : QWidget(parent) {
    visdata.numLayers = 0;
    visdata.minx = HUGE_VAL; visdata.miny = HUGE_VAL;
    visdata.maxx =-HUGE_VAL; visdata.maxy =-HUGE_VAL;
    antialiasing = true;
    mouseLeft = false;
    axisEqual = true;
    axiscolor = QColor("black");
    meshcolor = QColor(200,200,200);
    showColorbar = true;
}

MeshPlot::~MeshPlot() {
}

void
MeshPlot::clearData() {
    visdata.minx = HUGE_VAL; visdata.miny = HUGE_VAL;
    visdata.maxx =-HUGE_VAL; visdata.maxy =-HUGE_VAL;
    visdata.layers.clear();
    visdata.numLayers = 0;
}

void
MeshPlot::loadMesh(Mesh *mesh) {
    Layer *l = new Layer;
    l->loadMesh(mesh, meshcolor);
    visdata.layers.push_back(l);
    visdata.numLayers++;

    visdata.minx = std::min(visdata.minx, l->minx);
    visdata.miny = std::min(visdata.miny, l->miny);
    visdata.maxx = std::max(visdata.maxx, l->maxx);
    visdata.maxy = std::max(visdata.maxy, l->maxy);

    visdata.window_x0 = visdata.minx;
    visdata.window_y0 = visdata.miny;
    visdata.window_x1 = visdata.maxx;
    visdata.window_y1 = visdata.maxy;
    qDebug() << "loadMesh_zoom" << visdata.window_x0 << visdata.window_y0
              << visdata.window_x1 << visdata.window_y1;

}

void
MeshPlot::fixAspectRatio() {
    double ww = visdata.window_x1 - visdata.window_x0,
           wh = visdata.window_y1 - visdata.window_y0,
           xc = 0.5*(visdata.window_x1 + visdata.window_x0),
           yc = 0.5*(visdata.window_y1 + visdata.window_y0);

    double ARGui    = double(width()) / double(height()),
           ARWindow = ww / wh;

    if (ARGui > ARWindow) {
        visdata.window_x1 = xc + (ARGui / ARWindow) * 0.5 * ww;
        visdata.window_x0 = xc - (ARGui / ARWindow) * 0.5 * ww;
    } else {
        visdata.window_y1 = yc + (ARWindow / ARGui) * 0.5 * wh;
        visdata.window_y0 = yc - (ARWindow / ARGui) * 0.5 * wh;
    }
//    refreshPixmap();
}

void
MeshPlot::resizeEvent(QResizeEvent *) {
    // Match window aspect ratio to the aspect ratio of the window.
    if (axisEqual) {
        fixAspectRatio();
    }
    refreshPixmap();
}

// Translation of the window is started by pressing the left mouse
// button. At mousePressEvent, the window and the local coordinates
// of the mouse cursor are stored in mp_XXX variables and mouseLeft
// is set to true. Then at each mouseMoveEvent event with
// mouseLeft=true, the window is translated from the original window
// by an amount, which is proportional to the translation of the
// mouse cursor relative to its at position mousePressEvent.
// At mouseReleaseEvent, the variables mouseLeft is set to false.

void
MeshPlot::mousePressEvent(QMouseEvent *event) {
    if (event->button() == Qt::LeftButton) {
        QPoint p = event->pos();
        mp_x0 = p.x();
        mp_y0 = p.y();
        mp_wx0 = visdata.window_x0; mp_wy0 = visdata.window_y0;
        mp_wx1 = visdata.window_x1; mp_wy1 = visdata.window_y1;
        mouseLeft = true;
    }
}

void
MeshPlot::mouseReleaseEvent(QMouseEvent *event) {
    if (event->button() == Qt::LeftButton) {
        mouseLeft = false;
    }
}

void
MeshPlot::mouseMoveEvent(QMouseEvent *event) {
    if (mouseLeft) {
        QPoint p = event->pos();
        double x, y;
        coordLocalMesh(p.x(), p.y(), x, y);
        double ww = visdata.window_x1 - visdata.window_x0,
               wh = visdata.window_y1 - visdata.window_y0,
               dx = (p.x() - mp_x0)*ww/width(),
               dy = -(p.y() - mp_y0)*wh/height();

        visdata.window_x0 = mp_wx0 - dx;
        visdata.window_x1 = mp_wx1 - dx;
        visdata.window_y0 = mp_wy0 - dy;
        visdata.window_y1 = mp_wy1 - dy;
        refreshPixmap();
    }
}

void
MeshPlot::keyPressEvent(QKeyEvent *event) {
    double w = visdata.window_x1 - visdata.window_x0,
           h = visdata.window_y1 - visdata.window_y0;

    //    ug() << "keyPressEvent";
    switch (event->key()) {
    case Qt::Key_1:
        visdata.window_x0 = visdata.minx;
        visdata.window_x1 = visdata.maxx;
        visdata.window_y0 = visdata.miny;
        visdata.window_y1 = visdata.maxy;

        if (axisEqual) {
            fixAspectRatio();
        }
        refreshPixmap();
        break;
    case Qt::Key_Plus:
        zoom(0.1);
        refreshPixmap();
        break;
    case Qt::Key_Minus:
        zoom(-0.1);
        refreshPixmap();
        break;
    case Qt::Key_Left:
        visdata.window_x0 -= 0.1*w;
        visdata.window_x1 -= 0.1*w;
        refreshPixmap();
        break;
    case Qt::Key_Right:
        visdata.window_x0 += 0.1*w;
        visdata.window_x1 += 0.1*w;
        refreshPixmap();
        break;
    case Qt::Key_Down:
        visdata.window_y0 -= 0.1*h;
        visdata.window_y1 -= 0.1*h;
        refreshPixmap();
        break;
    case Qt::Key_Up:
        visdata.window_y0 += 0.1*h;
        visdata.window_y1 += 0.1*h;
        refreshPixmap();
        break;
    case Qt::Key_Escape:
//        exit(0);
        break;
    }
}

// This visualization code employs two coordinate systems:
// 1. The local coordinate system of the window with width()*height()
//    pixels.
// 2. The global coordinate system associated to the nodes of the mesh.
//    The following two functions are used for coordinate
//    transformations between these two systems.

void
MeshPlot::coordMeshLocal(double x, double y, double &lx, double &ly) {
    lx = width() * (x - visdata.window_x0)
                 / (visdata.window_x1 - visdata.window_x0);
    ly = height() * (1- (y - visdata.window_y0)
                 / (visdata.window_y1 - visdata.window_y0));
}

void
MeshPlot::coordLocalMesh(double lx, double ly, double &x, double &y) {
    x = visdata.window_x0
      + (lx / width())  * (visdata.window_x1 - visdata.window_x0);
    y = visdata.window_y0
      + (ly / height()) * (visdata.window_y1 - visdata.window_y0);
}

// The following three functions are wrappers for the drawing primitives
// in terms of global coordinates.

void
MeshPlot::drawPointGlobal(QPainter &painter, double x, double y, int size) {
    double lx, ly;
    coordMeshLocal(x, y, lx, ly);
    painter.drawRect((int)(lx-size/2), (int)(ly-size/2), size, size);
}

void
MeshPlot::drawLineGlobal(QPainter &painter, double x0, double y0,
                         double x1, double y1) {
    double lx0, ly0, lx1, ly1;
    coordMeshLocal(x0, y0, lx0, ly0);
    coordMeshLocal(x1, y1, lx1, ly1);

    painter.drawLine(lx0, ly0, lx1, ly1);
}

void
MeshPlot::drawPolyGlobal(QPainter &painter, QPolygonF P, QColor color) {
    int npoints = P.size();

    QPolygonF Pnew;
    Pnew.reserve(npoints);

    for (int ind_point = 0; ind_point < npoints; ind_point++) {
        QPointF p = P[ind_point];
        double lx, ly;
        coordMeshLocal(p.x(), p.y(), lx, ly);

        QPointF pnew(lx, ly);
        Pnew << pnew;
    }
    QPainterPath path;
    path.addPolygon(Pnew);
    painter.fillPath(path, QBrush(color));
//    painter.setBrush(QBrush(color));
//    painter.drawPolygon(Pnew);
}

// The mouse wheel is used to control the zoom. MeshPlot::zoom
// is called by MeshPlot::wheelEvent and scales the window size.

void
MeshPlot::zoom(double dsf) {
    double xc = 0.5*(visdata.window_x0 + visdata.window_x1),
           yc = 0.5*(visdata.window_y0 + visdata.window_y1),
           ww = visdata.window_x1 - visdata.window_x0,
           wh = visdata.window_y1 - visdata.window_y0;

    visdata.window_x0 = xc - 0.5*(1 - dsf)*ww;
    visdata.window_x1 = xc + 0.5*(1 - dsf)*ww;
    visdata.window_y0 = yc - 0.5*(1 - dsf)*wh;
    visdata.window_y1 = yc + 0.5*(1 - dsf)*wh;

    refreshPixmap();

}

void
MeshPlot::wheelEvent(QWheelEvent *event) {
//  QPoint p = event->pos();
    zoom(0.0005*event->delta());
}

// This implementation uses double buffering. That is, the graphical
// data is first drawn into a back buffer (QPixmap pixmap). At each
// paintEvent, the back buffer is then drawn into the front buffer
// using drawPixmap.

void
MeshPlot::paintEvent(QPaintEvent *) {
    QStylePainter painter(this);
    painter.drawPixmap(0, 0, pixmap);
}

void
MeshPlot::refreshPixmap() {
    pixmap = QPixmap(size());
    pixmap.fill();

    //qDebug() << visdata.window_x0 << visdata.window_y0 << visdata.window_x1 << visdata.window_y1;

    QPainter painter(&pixmap);
    painter.setRenderHint(QPainter::Antialiasing, antialiasing);
    if (axisEqual) {
        fixAspectRatio();
    }

    QPen penGrid, penMesh;
    penGrid.setColor(axiscolor);
    painter.setPen(penGrid);

    drawLineGlobal(painter, 0, visdata.miny, 0, visdata.maxy);
    drawLineGlobal(painter, visdata.minx, 0, visdata.maxx, 0);

    if (visdata.numLayers == 0) {
        painter.drawLine(0, height()/2, width(), height()/2);
        painter.drawLine(width()/2, 0, width()/2, height());
    }

    painter.setPen(penMesh);

    for (int ind_layer = 0; ind_layer < visdata.numLayers; ind_layer++) {
        Layer *layer = visdata.layers[ind_layer];

        if (layer->polysVisible) {
            for (int ind_poly = 0; ind_poly < layer->numPolygons; ind_poly++) {
                Polygon p = layer->polygons[ind_poly];
                drawPolyGlobal(painter, p.poly, p.color);
            }
        }
        if (layer->linesVisible) {
            painter.setPen(penMesh);
            //qDebug() << "nL" << layer->numLines;
            for (int ind_line = 0; ind_line < layer->numLines; ind_line++) {
                Line line = layer->lines[ind_line];
                penMesh.setColor(line.color);
                penMesh.setWidth(line.width);
                painter.setPen(penMesh);
                drawLineGlobal(painter, line.x0, line.y0, line.x1, line.y1);
            }
        }
        if (layer->nodesVisible) {
            for (int ind_node = 0; ind_node < layer->numNodes; ind_node++) {
                Node node = layer->nodes[ind_node];
                drawPointGlobal(painter, node.x, node.y, node.size);
            }
        }
    }
    if (showColorbar) {
        for (int ind_layer=0; ind_layer < visdata.numLayers; ind_layer++) {
            Layer *layer = visdata.layers[ind_layer];
            for (int ind_linegroup=0; ind_linegroup < layer->numLineGroups; ind_linegroup++) {
                QColor linecolor = layer->lineColormap[ind_linegroup];
                double dy = (double)height() / (double)(layer->numLineGroups + 1);
                int y = (int)(dy*(0.5+(double)ind_linegroup));

                penGrid.setColor(linecolor);
                penGrid.setWidth(2);
                painter.setPen(penGrid);

                double val = layer->lineValues[ind_linegroup];
                //qDebug() << val;
                QString vals;
                if (val == 0) {
                    vals.sprintf("0");
                } else {
                    int val_exp = (int)std::floor(std::log10(std::fabs(val)));
                    double val_coeff = val/std::pow(10, val_exp);
                    if (std::abs(val_exp) < 2) {
                        vals.sprintf("%.4f", val);
                    } else {
                        vals.sprintf("%.4fE%d", val_coeff, val_exp);
                    }
                }
                painter.fillRect(width()-20, y, 20, (int)dy+1, linecolor);
            //  painter.drawLine(width()-20, y, width(), y);
                penGrid.setColor(QColor("black"));
                painter.setPen(penGrid);

                painter.drawText(width()-130, y+4, 100, dy, Qt::AlignRight, vals);
            }
        }
    }
    update();
}

Layer::Layer(int lineReserve, int nodeReserve, int polygonReserve) {
    lines.reserve(lineReserve);
    nodes.reserve(nodeReserve);
    polygons.reserve(polygonReserve);
    numLines = 0;
    numNodes = 0;
    numPolygons = 0;
    numLineGroups = 0;

    nodesVisible = false;
    linesVisible = true;
    polysVisible = true;
    nodeLabelsVisible = false;
    lineLabelsVisible = false;
    polyLabelsVisible = false;
}

void
Layer::addNode(double x, double y, double size, QString label) {
    Node node;
    node.x = x;
    node.y = y;
    node.size = size;
    node.label = label;
    nodes.push_back(node);
    numNodes++;
}

void
Layer::addLine(double x0, double y0, double x1, double y1,
               double width, QString label, QColor color, double value) {
    Line line;
    line.x0 = x0;
    line.x1 = x1;
    line.y0 = y0;
    line.y1 = y1;
    line.width = width;
    line.label = label;
    line.color = color;
    line.value = value;
    lines.push_back(line);
    numLines++;
}

void
Layer::addLine(QPointF p0, QPointF p1, double width, QString label, QColor color, double value) {
    addLine(p0.x(), p0.y(), p1.x(), p1.y(), width, label, color, value);
}

void
Layer::addPolygon(QPolygonF poly, QString label, QColor color) {
    Polygon p;
    p.label = label;
    p.poly = poly;
    p.color = color;
    p.fillrule = Qt::OddEvenFill;

    polygons.push_back(p);
    numPolygons++;
}

void
Layer::reserveLines(int num) {
    lines.reserve(num);
}

void
Layer::reserveNodes(int num) {
    nodes.reserve(num);
}

void
Layer::reservePolygons(int num) {
    polygons.reserve(num);
}

// This code constructs a set of numContours contour lines from the
// degrees of freedom associated to the mesh nodes. DoFs associated
// to higher order basis functions are ignored in the visualization.
// Note: To avoid instability, the double *solution has to contain
// mesh->num_nodes values.

void
Layer::createContours(Mesh *mesh, QVector<double> &solution, int numContours) {
    maxvalue=-HUGE_VAL, minvalue=HUGE_VAL;
    Q_ASSERT(solution.size() == mesh->num_nodes);

    for (unsigned int ind=0; ind < mesh->num_nodes; ind++) {
        if (solution[ind] > maxvalue) {
            maxvalue = solution[ind];
        }
        if (solution[ind] < minvalue) {
            minvalue = solution[ind];
        }
    }
    //qDebug() << "min/max" << minvalue << maxvalue;
    double dval = (maxvalue - minvalue) / (numContours-1);

    // For each triangle and contour curve, we first find out the edges
    // of the triangle, which intersect with the curve. Then, a line
    // segment is added between the intersecting points.


    numLineGroups = numContours;
    for (unsigned int ind_triangle=0; ind_triangle < mesh->num_triangles; ind_triangle++) {
        Mesh_Element *elem = mesh->elements[mesh->triangles[ind_triangle]];
        // Indices of the triangle corner nodes.
        unsigned int node0 = elem->nodes[0]-1,
            node1 = elem->nodes[1]-1,
            node2 = elem->nodes[2]-1;

        // Values of the target function at the corner nodes
        double nodeval0 = solution[node0],
               nodeval1 = solution[node1],
               nodeval2 = solution[node2];
        // Locations of the target nodes in global coordinates.
        QPointF np0(mesh->nodes[node0*3], mesh->nodes[node0*3+1]),
                np1(mesh->nodes[node1*3], mesh->nodes[node1*3+1]),
                np2(mesh->nodes[node2*3], mesh->nodes[node2*3+1]);
        // Locations of the intersections.
        QPointF ne0, ne1, ne2;

        for (int ind_contour = 0; ind_contour < numContours; ind_contour++) {
//        for (int ind_contour = 10; ind_contour < 16; ind_contour++) {
            double cval = minvalue + ind_contour * dval;
            // Find out the intersecting nodes.
            bool edge0 = ((cval>nodeval0) && (cval<nodeval1)) ||
	      ((cval<nodeval0) && (cval>nodeval1));
            bool edge1 = ((cval>nodeval1) && (cval<nodeval2)) ||
	      ((cval<nodeval1) && (cval>nodeval2));
            bool edge2 = ((cval>nodeval2) && (cval<nodeval0)) ||
	      ((cval<nodeval2) && (cval>nodeval0));

            if (edge0) {
                ne0 = np0 + ((cval-nodeval0)/(nodeval1-nodeval0))*(np1-np0);
            }
            if (edge1) {
                ne1 = np1 + ((cval-nodeval1)/(nodeval2-nodeval1))*(np2-np1);
            }
            if (edge2) {
                ne2 = np2 + ((cval-nodeval2)/(nodeval0-nodeval2))*(np0-np2);
            }
            QColor lc;
            lc.setHsvF((cval-minvalue)/(maxvalue-minvalue), 1, 1);
            lineColormap[ind_contour] = lc;

            lineValues[ind_contour] = cval;

            if (edge0 && edge1) {
                addLine(ne0, ne1, 2, QString(), lc, cval);
                lineGroups[numLines] = ind_contour;
            }
            if (edge0 && edge2) {
                addLine(ne0, ne2, 2, QString(), lc, cval);
                lineGroups[numLines] = ind_contour;
            }
            if (edge1 && edge2) {
                addLine(ne1, ne2, 2, QString(), lc, cval);
                lineGroups[numLines] = ind_contour;
            }
        }
    }

}


// This code constructs a set of lines and nodes from a FE mesh.
void
Layer::loadMesh(Mesh *mesh, QColor meshcolor) {

    reserveLines(3*mesh->num_triangles);
    reserveNodes(mesh->num_nodes);

    double x, y;
    minx = HUGE_VAL; miny = HUGE_VAL;
    maxx =-HUGE_VAL; maxy =-HUGE_VAL;
    for (unsigned int ind_node = 0;ind_node < mesh->num_nodes; ind_node++) {
        x = mesh->nodes[ind_node*3];
        y = mesh->nodes[ind_node*3+1];
        addNode(x, y, 1, QString());

        //qDebug() << x << y;

        maxx = std::max(x, maxx);
        maxy = std::max(y, maxy);
        minx = std::min(x, minx);
        miny = std::min(y, miny);
    }
    //qDebug() << "Layer bnd:"<< minx << miny << maxx << maxy;

    for (unsigned int ind_triangle = 0; ind_triangle < mesh->num_triangles; ind_triangle++) {
        Mesh_Element *elem = mesh->elements[mesh->triangles[ind_triangle]];
        int node0 = elem->nodes[0]-1;
        int node1 = elem->nodes[1]-1;
        int node2 = elem->nodes[2]-1;

        double x0 = mesh->nodes[node0*3],
               y0 = mesh->nodes[node0*3+1],
               x1 = mesh->nodes[node1*3],
               y1 = mesh->nodes[node1*3+1],
               x2 = mesh->nodes[node2*3],
               y2 = mesh->nodes[node2*3+1];

        // Note that lines on interfaces between triangles are added twice.
        // Obviously, this is ugly and slow.
        addLine(x0, y0, x1, y1, 1, QString(), meshcolor);
        addLine(x1, y1, x2, y2, 1, QString(), meshcolor);
        addLine(x2, y2, x0, y0, 1, QString(), meshcolor);
    }
}

// This code colors physical regions with triangles of given color.
// The color associated to each physical region is given in the map
// physmap. When a physical region is not mapped to any color, the
// region is not colored.
void
Layer::colorPhysicals(Mesh *mesh, QMap<int, QColor> physmap) {
    for (unsigned int ind_triangle = 0; ind_triangle < mesh->num_triangles; ind_triangle++) {
        Mesh_Element *elem = mesh->elements[mesh->triangles[ind_triangle]];

        if (physmap.contains(elem->physical)) {
            int node0 = elem->nodes[0]-1;
            int node1 = elem->nodes[1]-1;
            int node2 = elem->nodes[2]-1;

            double x0 = mesh->nodes[node0*3],
                   y0 = mesh->nodes[node0*3+1],
                   x1 = mesh->nodes[node1*3],
                   y1 = mesh->nodes[node1*3+1],
                   x2 = mesh->nodes[node2*3],
                   y2 = mesh->nodes[node2*3+1];
            QPointF p0 = QPointF(x0, y0),
                    p1 = QPointF(x1, y1),
                    p2 = QPointF(x2, y2);
            QPolygonF p;
            p << p0 << p1 << p2;
            QColor color = physmap[elem->physical];
            addPolygon(p, QString(""), color);

        }
    }
}
