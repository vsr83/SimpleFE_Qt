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

#ifndef PHYSLIST_H
#define PHYSLIST_H

#include <QDialog>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QStack>
#include <QPushButton>
#include <QLabel>

#include <map>
#include "region.h"

class physList : public QDialog {
    Q_OBJECT
public:
    physList(std::map <int, Region> _physmap, QWidget *parent = 0);
    physList(QWidget *parent = 0);
    ~physList();

    void setFirst(std::map <int, Region> _physmapFirst);
    void updateList(std::map <int, Region> _physmap, bool toUndoStack);
signals:
    void physChanged(std::map<int, Region> _physmap, int phys=0);
    void physSelected(int phys);

    void resetClicked(std::map <int, Region> physmap);
    void undoClicked(std::map <int, Region> physmap);
    void redoClicked(std::map <int, Region> physmap);
private:
    QTableWidget *tableWidget;
    QHBoxLayout *layout;
    QVBoxLayout *toolbarLayout;
    QPushButton *undoButton, *redoButton, *clearButton, *resetButton;

    std::map <int, Region> physmap, physmapFirst;
    QStack <std::map <int, Region> > undoStack, redoStack;
    double mu0;
private slots:
    void cellChanged(int row, int col);
    void rowSelected(int row);
    void redo();
    void undo();
    void clearStacks();
    void reset();
};

#endif // PHYSLIST_H
