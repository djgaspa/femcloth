#include <iostream>
#include <vector>
#include <chrono>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <Eigen/IterativeLinearSolvers>
#include <QTimerEvent>
#include "Simulator.hpp"
extern "C" {
#define ANSI_DECLARATORS
#define REAL double
#define VOID void
#include "triangle.h"
}

static
float distance(float a, float b, float c, float x, float y)
{
    return a * x + b * y + c;
}

void Simulator::timerEvent(QTimerEvent* e)
{
    if (e->timerId() == timer_update) {
        QVector<float> ver(x.size());
        for (int i = 0; i < x.size(); ++i)
            ver[i] = x[i];
        emit updateGeometry(ver);
        return;
    }
    const int n_points = x.size() / 3;
    Eigen::Map<Eigen::MatrixXi> conn((int*)tri.data(), 3, tri.size() / 3);
    std::vector<Eigen::Triplet<float>> nzs, nzb, nz0;
    nzs.reserve(9 * 9 * conn.cols());
    nzb.reserve(9 * 9 * conn.cols());
    nz0.reserve(9 * 6 * conn.cols());
    for (int i = 0; i < conn.cols(); ++i) {
        const auto t = conn.col(i);
        Eigen::Vector3f v[3] {x.segment<3>(3 * t[0]), x.segment<3>(3 * t[1]), x.segment<3>(3 * t[2])};
        Eigen::Matrix3f Ri;
        Ri.col(0) = (v[0] - v[2]).normalized();
        const auto e2 = (v[1] - v[2]).normalized();
        Ri.col(1) = (e2 - Ri.col(0).dot(e2) * Ri.col(0)).normalized();
        Ri.col(2) = Ri.col(0).cross(Ri.col(1)).normalized();
        Eigen::Matrix<float, 3, 2> R = (Ri.transpose() * R0[i]).leftCols<2>();
        Eigen::Matrix<float, 2, 3> Rt = R.transpose();
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                const auto Q0 = R * Ks[a][b][i];
                const auto Qs = Q0 * Rt;
                //const auto Qb = (R * Kb[a][b][i]) * Rt; //TODO missing projections
                for (int k = 0; k < 3; ++k)
                    for (int l = 0; l < 3; ++l) {
                        nzs.push_back(Eigen::Triplet<float>(3 * t[a] + k, 3 * t[b] + l, Qs(k, l)));
                        //nzb.push_back(Eigen::Triplet<float>(3 * t[a] + k, 3 * t[b] + l, Qb(k, l)));
                        if (l < 2)
                            nz0.push_back(Eigen::Triplet<float>(3 * t[a] + k, 2 * t[b] + l, Q0(k, l)));
                    }
            }
    }
    Eigen::SparseMatrix<float> As(3 * n_points, 3 * n_points), Ab(3 * n_points, 3 * n_points), A0(3 * n_points, 2 * n_points);
    As.setFromTriplets(nzs.begin(), nzs.end());
    //Ab.setFromTriplets(nzb.begin(), nzb.end());
    A0.setFromTriplets(nz0.begin(), nz0.end());
    Eigen::SparseMatrix<float> B = 0.5f * As;
    std::vector<Eigen::Triplet<float>> nzm;
    nzm.reserve(3 * n_points);
    for (int i = 0; i < n_points; ++i) {
        const Eigen::Vector2f pos = x0.segment<2>(2 * i);
        for (int j = 0; j < 3; ++j)
            nzm.push_back(Eigen::Triplet<float>(3 * i + j, 3 * i + j, pos[1] > 0.9f ? 1e10f : mass[i]));
    }
    Eigen::SparseMatrix<float> M(3 * n_points, 3 * n_points);
    M.setFromTriplets(nzm.begin(), nzm.end());
    Eigen::VectorXf F = Eigen::VectorXf::Zero(3 * n_points);
    for (int i = 0; i < n_points; ++i) {
        const Eigen::Vector2f pos = x0.segment<2>(2 * i);
        F[3 * i + 1] = pos[1] > 0.9f ? 0.0f : -9.81f * mass[i];
        if (pos[0] > 0.4f && pos[0] < 0.6f && pos[1] < 0.1f) {
            //F[3 * i + 0] += -0.6f;
            //F[3 * i + 1] += 4.0f;
            //F[3 * i + 2] += 1.0f;
        }
    }
    const auto t1 = clock::now();
    const float dt = std::chrono::duration<float>(t1 - t).count();
    t = t1;
    Eigen::SparseMatrix<float> A = M + dt * dt * As + dt * B;
    Eigen::VectorXf b = M * v + dt * (F - As * x + A0 * x0);
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> solver;
    v = solver.compute(A).solveWithGuess(b, v);
    const auto info = solver.info();
    if (info == Eigen::ComputationInfo::NoConvergence) {
        std::cerr << "No convergence" << std::endl;
        return;
    }
    //std::cout << "solved with: " << solver.iterations() << " iterations" << std::endl;
    x += dt * v;
}

Simulator::Simulator(QObject *parent) :
    QObject(parent)
{
    qRegisterMetaType<QVector<unsigned>>("QVector<unsigned>");
    qRegisterMetaType<QVector<float>>("QVector<float>");
}

void Simulator::start()
{
    ::triangulateio in, out, vorout;
    std::vector<REAL> pointlist{0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    in.pointlist = pointlist.data();
    in.numberofpoints = pointlist.size() / 2;
    in.numberofpointattributes = 0;
    in.pointmarkerlist = nullptr;
    out.pointlist = nullptr;
    out.pointmarkerlist = nullptr;
    out.trianglelist = nullptr;
    out.neighborlist = nullptr;
    out.edgelist = nullptr;
    out.edgemarkerlist = nullptr;
    vorout.pointlist = nullptr;
    vorout.edgelist = nullptr;
    vorout.normlist = nullptr;
    std::string current_locale = std::setlocale(LC_ALL, nullptr);
    std::setlocale(LC_ALL, "C");
    ::triangulate(const_cast<char*>("QznevDq33a0.0018"), &in, &out, &vorout);
    std::setlocale(LC_ALL, current_locale.c_str());
    std::cout << "Number of points: " << out.numberofpoints << std::endl;
    std::vector<float> ver(3 * out.numberofpoints);
    for (int i = 0; i < out.numberofpoints; ++i) {
        ver[3 * i + 0] = out.pointlist[2 * i + 0];
        ver[3 * i + 1] = out.pointlist[2 * i + 1];
    }
    emit updateGeometry(QVector<float>::fromStdVector(ver));
    tri.resize(3 * out.numberoftriangles);
    std::copy(out.trianglelist, out.trianglelist + out.numberoftriangles * 3, tri.begin());
    emit updateConnectivity(QVector<unsigned>::fromStdVector(tri));
    std::vector<float> voronoi_ver(3 * vorout.numberofpoints);
    for (int i = 0; i < vorout.numberofpoints; ++i) {
        voronoi_ver[3 * i + 0] = vorout.pointlist[2 * i + 0];
        voronoi_ver[3 * i + 1] = vorout.pointlist[2 * i + 1];
    }
    std::vector<unsigned> voronoi_edges(2 * vorout.numberofedges);
    std::copy(vorout.edgelist, vorout.edgelist + 2 * vorout.numberofedges, voronoi_edges.begin());
    for (int i = 0; i < out.numberofedges; ++i) {
        if (out.edgemarkerlist[i] != 0) {
            const float x = vorout.pointlist[2 * vorout.edgelist[2 * i]];
            const float y = vorout.pointlist[2 * vorout.edgelist[2 * i] + 1];
            if (x > 1.0f || x < 0.0f || y < 0.0f || y > 1.0f) {
                voronoi_edges[2 * i + 1] = voronoi_edges[2 * i];
                continue;
            }
            const float nx = vorout.normlist[2 * i];
            const float ny = vorout.normlist[2 * i + 1];
            voronoi_ver.push_back(std::min(1.0f, std::max(0.0f, x + nx)));
            voronoi_ver.push_back(std::min(1.0f, std::max(0.0f, y + ny)));
            voronoi_ver.push_back(0);
            const int idx = voronoi_ver.size() / 3 - 1;
            voronoi_edges[2 * i + 1] = idx;
        }
        else {
            Eigen::Vector2f p0 {vorout.pointlist[2 * vorout.edgelist[2 * i]], vorout.pointlist[2 * vorout.edgelist[2 * i] + 1]};
            Eigen::Vector2f p1 {vorout.pointlist[2 * vorout.edgelist[2 * i + 1]], vorout.pointlist[2 * vorout.edgelist[2 * i + 1] + 1]};
            bool outside0 = p0[0] < 0.0f || p0[0] > 1.0f || p0[1] < 0.0f || p0[1] > 1.0f;
            bool outside1 = p1[0] < 0.0f || p1[0] > 1.0f || p1[1] < 0.0f || p1[1] > 1.0f;
            if (outside0 == false && outside1 == false)
                continue;
            unsigned& idx = voronoi_edges[2 * i + (outside0 == true ? 0 : 1)];
            const Eigen::Vector2f& q = outside0 == true ? p0 : p1;
            const float a = p1[1] - p0[1];
            const float b = p0[0] - p1[0];
            const float c = p0[0] * (p0[1] - p1[1]) + p0[1] * (p1[0] - p0[0]);
            float a1, b1, c1;
            if (distance(1, 0, 0, q[0], q[1]) < 0.0f) {
                a1 = 1; b1 = 0; c1 = 0;
            }
            else if (distance(1, 0, -1, q[0], q[1]) > 0.0f) {
                a1 = 1; b1 = 0; c1 = -1;
            }
            else if (distance(0, 1, 0, q[0], q[1]) < 0.0f) {
                a1 = 0; b1 = 1; c1 = 0;
            }
            else if (distance(0, 1, -1, q[0], q[1]) > 0.0f) {
                a1 = 0; b1 = 1; c1 = -1;
            }
            else {
                std::cerr << "Invalid connectivity" << std::endl;
                continue;
            }
            Eigen::Matrix2f A;
            A << a, b, a1, b1;
            const Eigen::Vector2f B(-c, -c1);
            const auto p = A.inverse() * B;
            voronoi_ver.push_back(p[0]);
            voronoi_ver.push_back(p[1]);
            voronoi_ver.push_back(0);
            idx = voronoi_ver.size() / 3 - 1;
        }
    }
    std::vector<std::list<std::pair<unsigned, unsigned>>> vt_list(out.numberofpoints);
    for (int i = 0; i < out.numberoftriangles; ++i) {
        const auto* t = &out.trianglelist[3 * i];
        for (unsigned j = 0; j < 3; ++j)
            vt_list[t[j]].push_back(std::make_pair(i, j));
    }
    std::vector<float> voronoi_areas(out.numberofpoints);
    for (int i = 0; i < out.numberofpoints; ++i) {
        Eigen::Vector2f point {out.pointlist[2 * i], out.pointlist[2 * i + 1]};
        std::vector<Eigen::Vector2f> region;
        auto& t_list = vt_list[i];
        if (out.pointmarkerlist[i] != 0) {
            region.push_back(point);
            auto it = t_list.begin();
            for (; it != t_list.end(); ++it) {
                const int t = it->first;
                const int e = it->second;
                const auto prev = out.neighborlist[3 * t + (e + 2) % 3];
                if (prev != -1)
                    continue;
                t_list.splice(t_list.begin(), t_list, it);
                const int e0 = out.trianglelist[3 * t + e], e1 = out.trianglelist[3 * t + (e + 1) % 3];
                int it_e = 0;
                for (; it_e < out.numberofedges; ++it_e)
                    if ((out.edgelist[2 * it_e] == e0 && out.edgelist[2 * it_e + 1] == e1) || (out.edgelist[2 * it_e] == e1 && out.edgelist[2 * it_e + 1] == e0))
                        break;
                if (it_e == out.numberofedges) {
                    std::cerr << "Invalid connectivity" << std::endl;
                    continue;
                }
                const int ve1 = voronoi_edges[2 * it_e + 1];
                const Eigen::Vector2f p {voronoi_ver[3 * ve1], voronoi_ver[3 * ve1 + 1]};
                if (p[0] >= 0.0f && p[0] <= 1.0f && p[1] >= 0.0f && p[1] <= 1.0f)
                    region.push_back(p);
                break;
            }
        }
        while (t_list.empty() == false) {
            const auto current_t = t_list.front().first;
            const auto current_e = t_list.front().second;
            Eigen::Vector2f p {vorout.pointlist[2 * current_t], vorout.pointlist[2 * current_t + 1]};
            if (p[0] < 0.0f || p[0] > 1.0f || p[1] < 0.0f || p[1] > 1.0f) {
                const auto next = out.trianglelist[3 * current_t + (current_e + 1) % 3];
                const auto prev = out.trianglelist[3 * current_t + (current_e + 2) % 3];
                if (out.pointmarkerlist[i] == 0) {
                    int it_e = 0;
                    for (; it_e < out.numberofedges; ++it_e)
                        if ((out.edgelist[2 * it_e] == i && out.edgelist[2 * it_e + 1] == next) || (out.edgelist[2 * it_e] == next && out.edgelist[2 * it_e + 1] == i))
                            break;
                    if (it_e == out.numberofedges) {
                        std::cerr << "Invalid connectivity" << std::endl;
                        continue;
                    }
                    int p_idx = voronoi_edges[2 * it_e + (vorout.edgelist[2 * it_e]  == current_t ? 0 : 1)];
                    region.push_back(Eigen::Vector2f{voronoi_ver[3 * p_idx], voronoi_ver[3 * p_idx + 1]});
                    for (it_e = 0; it_e < out.numberofedges; ++it_e)
                        if ((out.edgelist[2 * it_e] == i && out.edgelist[2 * it_e + 1] == prev) || (out.edgelist[2 * it_e] == prev && out.edgelist[2 * it_e + 1] == i))
                            break;
                    if (it_e == out.numberofedges) {
                        std::cerr << "Invalid connectivity" << std::endl;
                        continue;
                    }
                    p_idx = p_idx = voronoi_edges[2 * it_e + (vorout.edgelist[2 * it_e]  == current_t ? 0 : 1)];
                    region.push_back(Eigen::Vector2f{voronoi_ver[3 * p_idx], voronoi_ver[3 * p_idx + 1]});
                }
                else {
                    const auto next = out.neighborlist[3 * current_t + (current_e + 1) % 3];
                    const auto prev = out.neighborlist[3 * current_t + (current_e + 2) % 3];
                    if ((next == -1 && prev == -1) || (next != -1 && prev != -1)) {
                        std::cerr << "Invalid connectivity" << std::endl;
                        continue;
                    }
                    const auto e0 = next == -1 ? prev : next;
                    const auto e1 = current_t;
                    int it_e = 0;
                    for (; it_e < out.numberofedges; ++it_e)
                        if ((vorout.edgelist[2 * it_e] == e0 && vorout.edgelist[2 * it_e + 1] == e1) || (vorout.edgelist[2 * it_e] == e1 && vorout.edgelist[2 * it_e + 1] == e0))
                            break;
                    if (it_e == vorout.numberofedges) {
                        std::cerr << "Invalid connectivity" << std::endl;
                        continue;
                    }
                    int p_idx = voronoi_edges[2 * it_e + (vorout.edgelist[2 * it_e]  == current_t ? 0 : 1)];
                    region.push_back(Eigen::Vector2f{voronoi_ver[3 * p_idx], voronoi_ver[3 * p_idx + 1]});
                }
            }
            else
                region.push_back(p);
            const auto next = out.neighborlist[3 * current_t + (current_e + 1) % 3];
            t_list.erase(t_list.begin());
            if (next == -1) {
                const int t = current_t, e = (current_e + 2) % 3;
                const int e0 = out.trianglelist[3 * t + e], e1 = out.trianglelist[3 * t + (e + 1) % 3];
                int it_e = 0;
                for (; it_e < out.numberofedges; ++it_e)
                    if ((out.edgelist[2 * it_e] == e0 && out.edgelist[2 * it_e + 1] == e1) || (out.edgelist[2 * it_e] == e1 && out.edgelist[2 * it_e + 1] == e0))
                        break;
                if (it_e == out.numberofedges) {
                    std::cerr << "Invalid connectivity" << std::endl;
                    continue;
                }
                const int ve1 = voronoi_edges[2 * it_e + 1];
                const Eigen::Vector2f p {voronoi_ver[3 * ve1], voronoi_ver[3 * ve1 + 1]};
                if (p[0] >= 0.0f && p[0] <= 1.0f && p[1] >= 0.0f && p[1] <= 1.0f)
                    region.push_back(p);
                continue;
            }
            auto it = std::find_if(t_list.begin(), t_list.end(), [next](const std::pair<unsigned, unsigned>& pair) {
                return pair.first == next;
            });
            if (it == t_list.end())
                continue;
            t_list.splice(t_list.begin(), t_list, it);
        }
        auto area = 0.0f;
        for (unsigned i = 0; i < region.size(); ++i) {
            const int j = (i == region.size() - 1) ? 0 : i + 1;
            area += region[i][0] * region[j][1] - region[j][0] * region[i][1];
        }
        area /= 2.0f;
        voronoi_areas[i] = area;
    }
    mass.resize(out.numberofpoints);
    for (int i = 0; i < out.numberofpoints; ++i)
        mass[i] = voronoi_areas[i] * 0.23f;
    std::cout << "Total voronoi areas: " << std::accumulate(voronoi_areas.begin(), voronoi_areas.end(), 0.0f) << std::endl;
    assert(std::fabs(std::accumulate(voronoi_areas.begin(), voronoi_areas.end(), 0.0f) - 1) < 1e-4f);
    emit updateVoronoi(QVector<unsigned>::fromStdVector(voronoi_edges), QVector<float>::fromStdVector(voronoi_ver));
    x0 = Eigen::VectorXf(2 * out.numberofpoints);
    for (int i = 0; i < 2 * out.numberofpoints; ++i)
        x0[i] = out.pointlist[i];
    x = Eigen::VectorXf(3 * out.numberofpoints);
    for (int i = 0; i < ver.size(); ++i)
        x[i] = ver[i];
    v = Eigen::VectorXf::Zero(3 * out.numberofpoints);
    ::trifree(out.pointlist);
    ::trifree(out.pointmarkerlist);
    ::trifree(out.trianglelist);
    ::trifree(out.neighborlist);
    ::trifree(out.edgelist);
    ::trifree(out.edgemarkerlist);
    ::trifree(vorout.pointlist);
    ::trifree(vorout.edgelist);
    ::trifree(vorout.normlist);
    Eigen::Map<Eigen::MatrixXf> points(ver.data(), 3, ver.size() / 3);
    Eigen::Map<Eigen::MatrixXi> conn((int*)tri.data(), 3, tri.size() / 3);
    R0.resize(conn.cols());
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            Ks[i][j].resize(conn.cols());
            Kb[i][j].resize(conn.cols());
        }
    for (int i = 0; i < conn.cols(); ++i) {
        const auto t = conn.col(i);
        Eigen::Vector3f p[3] {points.col(t[0]), points.col(t[1]), points.col(t[2])};
        auto& R = R0[i];
        R.col(0) = (p[0] - p[2]).normalized();
        Eigen::Vector3f e2 = (p[1] - p[2]).normalized();
        R.col(1) = (e2 - R.col(0).dot(e2) * R.col(0)).normalized();
        R.col(2) = R.col(0).cross(R.col(1));
        Eigen::Matrix3f m;
        m << 1, p[0].segment<2>(0).transpose(), 1, p[1].segment<2>(0).transpose(), 1, p[2].segment<2>(0).transpose();
        const auto det = m.determinant();
        const auto area = det / 2.0;
        Eigen::Matrix<float, 2, 3> N; // form factor
        N << p[1].y() - p[2].y(), p[2].y() - p[0].y(), p[0].y() - p[1].y(), p[2].x() - p[1].x(), p[0].x() - p[2].x(), p[1].x() - p[0].x();
        N /= det;
        std::array<float, 16> C;
        C.fill(0.0f);
        C[0] = 245;
        C[15] = 366;
        C[5] = C[10] = C[9] = C[6] = 0.38;
        C[3] = C[12] = 61.1;
        std::array<float, 2> B {0.013, 0.037};
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                Eigen::Matrix2f& K = Ks[a][b][i];
                K.setZero();
                for (int i = 0; i < 2; ++i)
                    for (int j = 0; j < 2; ++j) {
                        float z = 0;
                        for (int k = 0; k < 2; ++k)
                            for (int l = 0; l < 2; ++l) {
                                int idx = l + 2 * k + 4 * j + 8 * i;
                                z += N(k, a) * C[idx] * N(l, b);
                            }
                        K(i, j) = z;
                    }
                K *= area;
            }
        for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b) {
                Eigen::Matrix2f& K = Kb[a][b][i];
                K << N(0, a) * B[0] * N(0, b), 0, 0, N(1, a) * B[1] * N(1, b);
                K *= area;
            }
    }
    timer_simulate = startTimer(0);
    timer_update = startTimer(33);
    t = clock::now();
}
