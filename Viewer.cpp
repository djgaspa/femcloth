#include <QThread>
#include <GL/glew.h>
#include "Viewer.hpp"
#include "Simulator.hpp"

extern "C" {
#define ANSI_DECLARATORS
#define REAL double
#define VOID void
#include "triangle.h"
}

Viewer::Viewer(QWidget *parent) :
    QGLViewer(parent),
    t(new QThread(this))
{
}

Viewer::~Viewer()
{
    t->quit();
    t->wait();
    ::glDeleteBuffers(vbo.size(), vbo.data());
    ::glDeleteVertexArrays(vao.size(), vao.data());
}

void Viewer::init()
{
    ::glewInit();
    ::glGenVertexArrays(vao.size(), vao.data());
    ::glGenBuffers(vbo.size(), vbo.data());
    ::glBindVertexArray(vao[0]);
    ::glEnableClientState(GL_VERTEX_ARRAY);
    ::glBindVertexArray(0);
    auto* s = new Simulator;
    s->moveToThread(t);
    connect(t, SIGNAL(finished()), s, SLOT(deleteLater()));
    connect(s, SIGNAL(updateConnectivity(QVector<uint>)), this, SLOT(updateConnectivity(QVector<uint>)));
    connect(s, SIGNAL(updateGeometry(QVector<float>)), this, SLOT(updateGeometry(QVector<float>)));
    connect(s, SIGNAL(updateVoronoi(QVector<uint>,QVector<float>)), this, SLOT(updateVoronoi(QVector<uint>,QVector<float>)));
    connect(t, SIGNAL(started()), s, SLOT(start()));
    t->start();
}

void Viewer::draw()
{
    ::glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    ::glColor3f(1.0f, 0.0f, 0.0f);
    ::glBindVertexArray(vao[0]);
    ::glDrawElements(GL_TRIANGLES, n_elements, GL_UNSIGNED_INT, 0);
    ::glBindVertexArray(vao[1]);
    ::glColor3f(0.0f, 1.0f, 0.0f);
    ::glDrawElements(GL_LINES, n_voronoi_elements, GL_UNSIGNED_INT, 0);
    ::glBindVertexArray(0);
    ::glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void Viewer::updateConnectivity(QVector<unsigned> tri)
{
    ::glBindVertexArray(vao[0]);
    ::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[0]);
    ::glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned) * tri.size(), tri.data(), GL_STATIC_DRAW);
    ::glBindVertexArray(0);
    n_elements = tri.size();
}

void Viewer::updateGeometry(QVector<float> ver)
{
    ::glBindBuffer(GL_ARRAY_BUFFER, vbo[1]);
    ::glBufferData(GL_ARRAY_BUFFER, sizeof(float) * ver.size(), ver.data(), GL_STATIC_DRAW);
    ::glBindVertexArray(vao[0]);
    ::glVertexPointer(3, GL_FLOAT, 0, 0);
    ::glBindVertexArray(0);
    ::glBindBuffer(GL_ARRAY_BUFFER, 0);
    updateGL();
}

void Viewer::updateVoronoi(const QVector<unsigned> edges, const QVector<float> ver)
{
    ::glBindVertexArray(vao[1]);
    ::glEnableClientState(GL_VERTEX_ARRAY);
    ::glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo[2]);
    ::glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned) * edges.size(), edges.data(), GL_STATIC_DRAW);
    ::glBindBuffer(GL_ARRAY_BUFFER, vbo[3]);
    ::glBufferData(GL_ARRAY_BUFFER, sizeof(float) * ver.size(), ver.data(), GL_STATIC_DRAW);
    ::glVertexPointer(3, GL_FLOAT, 0, 0);
    ::glBindBuffer(GL_ARRAY_BUFFER, 0);
    ::glBindVertexArray(0);
    n_voronoi_elements = edges.size();
}
