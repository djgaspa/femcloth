#pragma once
#include <array>
#include <vector>
#include <QGLViewer/qglviewer.h>
#include <Eigen/Core>
#include <Eigen/StdVector>

class QTimerEvent;

class Viewer : public QGLViewer
{
    Q_OBJECT

    std::array<unsigned, 2> vao;
    std::array<unsigned, 4> vbo;
    unsigned n_elements = 0, n_voronoi_elements = 0;
    std::vector<float> ver;
    std::vector<int> tri;
    std::vector<float> mass;
    std::array<std::array<std::vector<Eigen::Matrix2f, Eigen::aligned_allocator<Eigen::Matrix2f>>, 3>, 3> Ks, Kb;
    std::vector<Eigen::Matrix3f, Eigen::aligned_allocator<Eigen::Matrix3f>> R0;

    Eigen::VectorXf x0, x, v;

    void timerEvent(QTimerEvent*) override;
public:
    explicit Viewer(QWidget *parent = 0);
    ~Viewer();
    
signals:
    
public slots:
    virtual void init() override;
    virtual void draw() override;
};
