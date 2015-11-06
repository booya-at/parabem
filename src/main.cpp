#include <iostream>
#include "vector.h"
#include "panel3.h"
#include "element_influence.h"
#include "case3.h"

    using namespace std;

    int main(int argc, char **argv) {
    PanelVector3 a(1,0, 1);
    PanelVector3 b(0,0,1);
    PanelVector3 c(-1,0,0);
    PanelVector3 d(0,3,4);
    Panel3 pan;
    pan.append_point(&a);
    pan.append_point(&b);
    pan.append_point(&c);
    pan.append_point(&d);
    pan.calc_geo();
    SymmetricPanel3* g = new SymmetricPanel3(&pan, a, b);
    g->calc_geo();
}