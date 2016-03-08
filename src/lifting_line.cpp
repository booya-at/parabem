#include "lifting_line.h"


LineSegment::LineSegment(Vector3& v1, Vector3& v2, Vector3& v_inf)
{
    this->v1 = &v1;
    this->v2 = &v2;
    this->mid = (v1 + v2) / 2;
    this->n = v_inf.cross(*this->v2 - *this->v1);
    this->n.normalize();
    this->t = v_inf;
    this->t.normalize();
    this->lift_factor = (*this->v2 - *this->v1).cross(Vector3(0,0,1)).norm();
}

double LineSegment::ind_influence(LineSegment& ls, double gamma)
{
    return gamma * (this->ind_influence(ls.mid, gamma).dot(this->n));
}

Vector3 LineSegment::ind_influence(Vector3& point, double gamma)
{
    return (vortex_3_0_half_infinity_v(*this->v2, this->t, point) -
            vortex_3_0_half_infinity_v(*this->v1, this->t, point)
            ) * gamma;
}

double LineSegment::ca()
{
    return this->best_gamma * this->b();
}

double LineSegment::b()
{
    return (*this->v1 - *this->v2).norm();
}

void LiftingLine::append_point(Vector3& v)
{
    this->lifting_line.push_back(v);
}


void LiftingLine::initialize(Vector3& v_inf)
{
    this->v_inf = v_inf;
    for (int i = 0; i < this->lifting_line.size() - 1; i++){
        LineSegment *ll = new LineSegment(this->lifting_line[i], this->lifting_line[i+1], this->v_inf);
        this->segments.push_back(ll);
    }
}

void LiftingLine::best_gamma(double cL)
{
    int matsize = this->segments.size() + 1;
    Eigen::MatrixXd matrix;
    Eigen::VectorXd rhs;
    Eigen::VectorXd result;
    matrix.resize(matsize, matsize);
    matrix.setZero();
    rhs.resize(matsize);
    rhs.setZero();
    result.resize(matsize);
    result.setZero();
    //   adding the main min-function
    for (int i = 0; i < (matsize - 1); i++)
    {
        for (int j = 0; j < (matsize - 1); j++)
        {
            matrix(i, j) += this->segments[i]->b() * this->segments[j]->ind_influence(*this->segments[i]);
            matrix(i, j) += this->segments[j]->b() * this->segments[i]->ind_influence(*this->segments[j]);
        }
    //     adding lagrange multiplicator
        matrix(i, matsize - 1) += this->segments[i]->lift_factor;
    }
    for (int i = 0; i < (matsize -1); i++)
    {
        matrix(matsize - 1, i) += this->segments[i]->lift_factor;
    }
    rhs(matsize - 1) += cL;
    
    result = matrix.fullPivHouseholderQr().solve(rhs);
    for (int i = 0; i < matsize - 1; i++)
    {
        this->segments[i]->best_gamma = result[i];
    }
}
