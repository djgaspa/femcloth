#pragma once
#include <vector>
#include <array>
#include <QObject>
#include <Eigen/Core>
#include <Eigen/StdVector>
#include <QVector>
#include <QMetaType>

class Simulator : public QObject
{
    Q_OBJECT

    std::vector<unsigned> tri;
    std::vector<float> mass;
    std::array<std::array<std::vector<Eigen::Matrix2f, Eigen::aligned_allocator<Eigen::Matrix2f>>, 3>, 3> Ks, Kb;
    std::vector<Eigen::Matrix3f, Eigen::aligned_allocator<Eigen::Matrix3f>> R0;

    Eigen::VectorXf x0, x, v;

    int timer_update, timer_simulate;

    void timerEvent(QTimerEvent*) override;

public:
    explicit Simulator(QObject *parent = 0);

signals:
    void updateConnectivity(const QVector<unsigned> tri);
    void updateGeometry(const QVector<float> ver);
    void updateVoronoi(const QVector<unsigned> lines, const QVector<float> ver);

public slots:
    void start();

};

Q_DECLARE_METATYPE(QVector<unsigned>)
Q_DECLARE_METATYPE(QVector<float>)
