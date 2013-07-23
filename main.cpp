#include <iostream>
#include <vector>
#include <QApplication>
#include "Viewer.hpp"

int main(int argc, char** argv)
{
    QApplication app(argc, argv);
    Viewer v;
    v.show();
    return app.exec();
}
