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

// This is an implementation of a simple Qt widget, which can be used
// for visualization of Finite Element (FE) meshes and scalar data.
//
// Graphical data is assembled into a vector of layers, where each
// layer is assembled from vectors of Node, Line and Polygon objects.
// Coordinates of all objects are expressed in global coordinates.
// The VisData object contains the vector of layers and a window,
// which determines the rectangular area drawn when the widget is
// painted. The window coordinates are changed by operations used for
// translation and zooming.

#ifndef MESHPLOT_H
#define MESHPLOT_H

#include <QWidget>
#include <QPen>
#include <QPixmap>
#include <QColor>
#include <QMap>
#include "mesh.h"

class Polygon;
class Line;
class Node;
class Layer;
class VisData;

class VisData {
public:
    QVector <Layer *> layers;
    int numLayers;

    double minx, miny, maxx, maxy;
    double window_x0, window_x1, window_y0, window_y1;
};

class MeshPlot : public QWidget
{
    Q_OBJECT

public:
    MeshPlot(QWidget *parent = 0);
    ~MeshPlot();
    void refreshPixmap();
    void fixAspectRatio();

    void loadMesh(Mesh *filename);
    void setAntiAliasing(bool aa) {antialiasing = aa;}
    void setAxisEqual   (bool ae) {axisEqual = ae;}
    void setShowColorbar(bool sc) {showColorbar = sc;}
    void clearData();
    void setAxisColor   (QColor c){axiscolor = c;}
    void setMeshColor   (QColor c){meshcolor = c;}

    VisData visdata;
public slots:
protected:
    void paintEvent       (QPaintEvent *event);
    void resizeEvent      (QResizeEvent *event);
    void mousePressEvent  (QMouseEvent *event);
    void mouseMoveEvent   (QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void keyPressEvent    (QKeyEvent *event);
    void wheelEvent       (QWheelEvent *event);
private: // Note: QWidget inherits width() and height().
    void drawPointGlobal  (QPainter &painter, double x, double y, int size);
    void drawLineGlobal   (QPainter &painter, double x0, double y0,
                           double x1, double y1);
    void drawPolyGlobal   (QPainter &painter, QPolygonF p, QColor color);
    void coordLocalMesh   (double lx, double ly, double &x, double &y);
    void coordMeshLocal   (double x, double y, double &lx, double &ly);

    QPixmap pixmap;      // The pixmap used as the back buffer.
    bool antialiasing;   // Use anti-aliasing.
    bool mouseLeft;      // Is left mouse button down?
    bool axisEqual;
    bool showColorbar;

    void zoom(double dsf);

    QColor meshcolor, axiscolor;

    // Temporary variables used to store the window and the mouse
    // location when left button is pressed to start translation of
    // the window.
    int mp_x0, mp_y0;
    double mp_wx0, mp_wx1, mp_wy0, mp_wy1;
};

class Line {
public:
    double x0, y0, x1, y1;
    double width;
    double value;
    QString label;
    QColor color;
};

class Node {
public:
    double x, y;
    double size;
    QString label;
};

class Polygon {
public:
    QPolygonF poly;
    QColor color;
    QString label;
    Qt::FillRule fillrule;
};

class Layer {
public:
    Layer(int lineReserve=0, int nodeReserve=0, int polygonReserve=0);

    void reserveLines   (int num);
    void reserveNodes   (int num);
    void reservePolygons(int num);

    // Functions for adding primitimes to layers.
    void addLine   (double x0, double y0, double x1, double y1, double width,
                    QString label, QColor color, double value=0);
    void addLine   (QPointF p0, QPointF p1, double width, QString label,
                    QColor color, double value=0);
    void addNode   (double x, double y, double size, QString label);
    void addPolygon(QPolygonF poly, QString label, QColor color);

    // Functions for creating visualization data.
    void loadMesh      (Mesh *mesh, QColor meshcolor);
    void colorPhysicals(Mesh *mesh, QMap<int, QColor>physmap);
    void createContours(Mesh *mesh, QVector<double> &solution, int numContours);

    QVector<Line> lines;
    QVector<Node> nodes;
    QVector<Polygon> polygons;

    QMap <int, int>    lineGroups;
    QMap <int, QColor> lineColormap;
    QMap <int, double> lineValues;
    int    numLineGroups;
    double minvalue, maxvalue;

    QString name;
    int numLines, numNodes, numPolygons;
    bool nodesVisible, linesVisible, polysVisible,
         nodeLabelsVisible, lineLabelsVisible, polyLabelsVisible;

    double minx, maxx, miny, maxy;

    QPen lineStyle;
};

#endif // MESHPLOT_H
