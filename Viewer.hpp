#pragma once
#include <array>
#include <vector>
#include <QGLViewer/qglviewer.h>
#include <QVector>

class QThread;

class Viewer : public QGLViewer
{
    Q_OBJECT

    QThread* t;

    std::array<unsigned, 2> vao;
    std::array<unsigned, 4> vbo;
    unsigned n_elements = 0, n_voronoi_elements = 0;

public:
    explicit Viewer(QWidget *parent = 0);
    ~Viewer();
    
signals:
    
public slots:
    virtual void init() override;
    virtual void draw() override;

private slots:
    void updateConnectivity(QVector<unsigned> tri);
    void updateGeometry(QVector<float> ver);
    void updateVoronoi(const QVector<unsigned> edges, const QVector<float> ver);
};
