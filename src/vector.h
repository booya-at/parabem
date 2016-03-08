#ifndef vector_H
#define vector_H

#include <Eigen/Core>
#include <vector>
#include <iostream>
#include <sstream>


//forward decleration
class Panel3;
class Panel2;


/**
 * A 2D Vector from the Eigen library -> http://eigen.tuxfamily.org/index.php?title=Main_Page.
*/

typedef Eigen::Vector2d Vector2;


/**
 * A 3D Vector from the Eigen library -> http://eigen.tuxfamily.org/index.php?title=Main_Page.
*/

typedef Eigen::Vector3d Vector3;


/** 
 * A Panel Vector2 is a Vector2 with some other members like a scalar potential, an vector 
 * velocity and a list of Panels which contain all the panels which use this vertex.
*/

class PanelVector2: public Vector2{
public:
    PanelVector2():Vector2(){}
    PanelVector2(double x, double y):Vector2(x, y){}
    template<typename OtherDerived> PanelVector2(const Eigen::MatrixBase<OtherDerived>& other):Vector2(other){};
    template<typename OtherDerived> PanelVector2& operator=(const Eigen::MatrixBase <OtherDerived>& other){
        this->Vector2::operator=(other);
        return *this;
    }
    int nr = -1;
    double potential = 0.;
    bool wake_vertex = false;
    Vector2 velocity = Vector2(0, 0);
    std::vector<Panel2*> panels;
    void reset_properties(){this->potential=0; this->velocity.setZero();};
};


/** 
 * A Panel Vector3 is a Vector3 with some other members like a scalar potential, an vector 
 * velocity and a list of Panels which contain all the panels which use this vertex.
*/

class PanelVector3: public Vector3{
public:
    PanelVector3():Vector3(){}
    PanelVector3(double x, double y, double z):Vector3(x, y, z){}
    template<typename OtherDerived> PanelVector3(const Eigen::MatrixBase<OtherDerived>& other):Vector3(other){};
    template<typename OtherDerived> PanelVector3& operator=(const Eigen::MatrixBase <OtherDerived>& other){
        this->Vector3::operator=(other);
        return *this;
    }
    
    int owner = 0;
    int nr = -1;
    double potential = 0.;
    double vorticity = 0;
    double cp = 0.;
    Vector3 velocity = Vector3(0, 0, 0);
    bool wake_vertex = false;
    std::vector<Panel3*> panels;		//all Panels containing this vertex
    void reset_properties(){this->potential=0; this->vorticity=0; this->velocity.setZero();this->cp=0;};
};

static Vector2 normal2(Vector2 direction)
{
    return Vector2(-direction.y(), direction.x());
};


#endif // vector_H
