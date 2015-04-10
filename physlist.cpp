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

#include "physlist.h"
#include <map>
#include <QHeaderView>
#include <QTableWidgetItem>
#include <QDebug>
#include <math.h>

void
physList::clearStacks() {
    qDebug() << "clearStacks()";
    undoStack.clear();
    redoStack.clear();
    updateList(physmap, false);
}

void
physList::setFirst(std::map<int, Region> _physmapFirst) {
    physmapFirst = _physmapFirst;
    reset();
}

void
physList::reset() {
    qDebug() << "reset()";
    undoStack.clear();
    redoStack.clear();
    updateList(physmapFirst, false);
    emit resetClicked(physmapFirst);
}

void
physList::redo() {
    qDebug() << "redo()";
    undoStack.push(physmap);
    std::map <int, Region> tmpPhys = redoStack.pop();
    updateList(tmpPhys, false);
    emit redoClicked(tmpPhys);
}

void
physList::undo() {
    qDebug() << "undo()";
    redoStack.push(physmap);
    std::map <int, Region> tmpPhys = undoStack.pop();
    updateList(tmpPhys, false);
    emit undoClicked(tmpPhys);
}

void
physList::rowSelected(int row) {
    QString physstr = tableWidget->item(row, 3)->text();
    bool ok;
    int phys = physstr.toInt(&ok);
    Q_ASSERT(ok);
    emit physSelected(phys);
    qDebug() << "row" << row << "selected";
}

void
physList::cellChanged(int row, int col) {
    qDebug() << "foo" << row << col << tableWidget->rowCount() << tableWidget->columnCount();
    qDebug() << tableWidget->item(row, 3);

    if (!tableWidget->item(row, 3)) return;

    QString celltext = tableWidget->item(row, col)->text();
    QString physstr  = tableWidget->item(row, 3)->text();
    qDebug() << "cellChanged" << row << col << celltext;

    bool ok;
    int phys = physstr.toInt(&ok);
    if (!ok) return;

    double cellvalue = celltext.toDouble(&ok);

    QString newtext;
    tableWidget->blockSignals(true);
    switch (col) {
    case 0:
        undoStack.push(physmap);
        redoStack.clear();

        physmap[phys].name = std::string(qPrintable(celltext));
        break;
    case 1:
        if (ok) {
            undoStack.push(physmap);
            redoStack.clear();
            qDebug() << "Stack1" << undoStack.size();

            physmap[phys].nu = 1/(mu0*cellvalue);
            emit physChanged(physmap, phys);
            newtext.sprintf("%.1f", cellvalue);
        } else {
            qDebug() << "invalid value" << celltext<< "-> revert back to" << physmap[phys].nu;
            newtext.sprintf("%.1f", 1/(mu0*physmap[phys].nu));
        }
        tableWidget->item(row, col)->setText(newtext);
        break;
    case 2:
        if (ok) {
            undoStack.push_back(physmap);
            redoStack.clear();
            qDebug() << "Stack2" << undoStack.size();

            physmap[phys].Js = cellvalue*1e6;
            emit physChanged(physmap, phys);
            newtext.sprintf("%.1f", cellvalue);
        } else {
            qDebug() << "invalid value" << physmap[phys].Js;
            newtext.sprintf("%.1f", physmap[phys].Js/1e6);
        }
        tableWidget->item(row, col)->setText(newtext);
        break;
    }
    tableWidget->blockSignals(false);

    QString undoString, redoString;
    undoString = undoString.sprintf("Undo (%d)", undoStack.size());
    redoString = redoString.sprintf("Redo (%d)", redoStack.size());
    undoButton->setText(undoString);
    redoButton->setText(redoString);

    undoButton->setEnabled(undoStack.size() > 0);
    redoButton->setEnabled(redoStack.size() > 0);
    clearButton->setEnabled((redoStack.size() > 0) || (undoStack.size() > 0));
}

void
physList::updateList(std::map<int, Region> _physmap, bool toUndoStack) {
    qDebug() << "update";
    if (toUndoStack) {
        undoStack.push_back(physmap);
    }
    physmap = _physmap;

    for (int ind_row=tableWidget->rowCount()-1; ind_row >= 0; ind_row--) {
        tableWidget->removeRow(ind_row);
    }

    int row = 0;
    for (std::map<int, Region>::iterator it = physmap.begin();
         it != physmap.end(); ++it) {
        int physical  = it->first;
        Region region = it->second;

        qDebug() << physical << region.name.data() << region.nu << region.Js;
        QTableWidgetItem *item_phys = new QTableWidgetItem,
                         *item_name = new QTableWidgetItem,
                         *item_nu   = new QTableWidgetItem,
                         *item_Js   = new QTableWidgetItem;
        tableWidget->insertRow(row);

        tableWidget->setItem(row, 0, item_name);
        tableWidget->setItem(row, 1, item_nu);
        tableWidget->setItem(row, 2, item_Js);
        tableWidget->setItem(row, 3, item_phys);
        item_phys->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);

        QString string_phys, string_label, string_nu, string_Js;
        string_phys = string_phys.sprintf("%d", physical);
        string_label= QString(region.name.data());
        string_nu   = string_nu.sprintf("%.1f", 1/(region.nu*mu0));
        string_Js   = string_Js.sprintf("%.1f", region.Js/1e6);

        item_phys->setFlags(item_phys->flags() ^ Qt::ItemIsEditable);

        tableWidget->item(row, 0)->setText(string_label);
        tableWidget->item(row, 1)->setText(string_nu);
        tableWidget->item(row, 2)->setText(string_Js);
        tableWidget->item(row, 3)->setText(string_phys);
        row++;
    }
    setWindowTitle(tr("Material Parameters of Physical Regions"));
}

physList::physList(std::map<int, Region> _physmap, QWidget *parent) : QDialog(parent) {
    physmapFirst = physmap = _physmap;
    tableWidget = new QTableWidget(0, 4);

    mu0 = 4*M_PI*1e-7;
    resize(600, 300);
    updateList(physmap, false);

    QString mustring = QString::fromUtf8("\u03BC (rel.)");
    QString Jsstring = QString::fromUtf8("J\u209B (A/mm\u00B2)");
    tableWidget->setHorizontalHeaderLabels(QStringList()
                                           << tr("Name")
                                           << mustring
                                           << Jsstring
                                           << tr("Physical"));
    tableWidget->horizontalHeader()->setSectionResizeMode(0, QHeaderView::Stretch);

    layout = new QHBoxLayout;
    layout->addWidget(tableWidget);
    setLayout(layout);

    toolbarLayout = new QVBoxLayout;
    layout->addLayout(toolbarLayout);
    undoButton  = new QPushButton(tr("Undo"));
    redoButton  = new QPushButton(tr("Redo"));
    clearButton = new QPushButton(tr("Clear"));
    resetButton = new QPushButton(tr("Reset"));

    redoButton->setEnabled(false);
    undoButton->setEnabled(false);
    clearButton->setEnabled(false);
    resetButton->setEnabled(true);

    connect(redoButton,  SIGNAL(clicked()), this, SLOT(redo()));
    connect(undoButton,  SIGNAL(clicked()), this, SLOT(undo()));
    connect(clearButton, SIGNAL(clicked()), this, SLOT(clearStacks()));
    connect(resetButton, SIGNAL(clicked()), this, SLOT(reset()));

    toolbarLayout->addWidget(undoButton);
    toolbarLayout->addWidget(redoButton);
    toolbarLayout->addWidget(clearButton);
    toolbarLayout->addWidget(resetButton);
    toolbarLayout->addStretch();

    connect(tableWidget, SIGNAL(cellChanged(int,int)), this, SLOT(cellChanged(int,int)));

    tableWidget->horizontalHeader()->setSectionsClickable(false);
    connect(tableWidget->verticalHeader(), SIGNAL(sectionClicked(int)), this, SLOT(rowSelected(int)));

    updateList(physmap, false);
}

physList::physList(QWidget *parent) {
    physList(physmap, parent);
}

physList::~physList() {

}
