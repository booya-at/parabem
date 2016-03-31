#ifndef BOOYA_JSON
#define BOOYA_JSON


#include <jsoncpp/json/json.h>
#include <vector>
#include "panel3.h"
#include "vector.h"
#include "case3.h"

// not ready yet!!!


using namespace std;
void case_from_json(ifstream &inputfile, Case3* c){
    Json::Value root;
    Json::Reader reader;
    vector<Panel3*> panels;

    bool parsingSuccessful = reader.parse( inputfile, root );
    if ( !parsingSuccessful )
    {
        // report to the user the failure and their locations in the document.
        cout  << "Failed to parse inputfile\n";
        return 1;
    }
    std::string encoding = root.get("encoding", "UTF-8" ).asString();

    const Json::Value json_config = root["config"];
    const Json::Value json_nodes = root["nodes"];
    const Json::Value json_panels = root["panels"];

    if (json_nodes.empty() || json_panels.empty()){
            cout << "Please Provide a 'panels', 'nodes' and 'config' section in json file" << endl;
            cout << keywords << endl;
            return 1;
    }

    config->node_number = json_nodes.size();
    config->panel_number = json_panels.size();
    cout << "Reading:\n\t" << config->node_number << " Nodes\n\t" << config->panel_number << " Panels" << endl;

    //INITIALISATION
    *nodes = (Eigen::Vector3d*)malloc(sizeof(Eigen::Vector3d) * config->node_number);
    *panels = (Panel*)malloc(sizeof(Panel) * config->panel_number);

    //V_INF
Json::Value json_v_inf = json_config["v_inf"];
config->v_inf = Eigen::Vector3d();
for (int index=0; index<3; index++)
    config->v_inf[index] = json_v_inf[index].asDouble();

    //READ NODES
    for (int i=0; i<config->node_number; i++) {
    (*nodes)[i] = Eigen::Vector3d();
    for (int index = 0; index < 3; index++) {
        (*nodes)[i][index] = json_nodes[i][index].asDouble();
    }
}
    //READ PANELS
    int j=0;
    for (int i=0; i<config->panel_number; i++){
            Json::Value json_thispanel = json_panels[i];
            Panel *thispanel = *panels + i;
    *thispanel = Panel();
            //cout << i << endl;
            // make assert (->dict)
            //NODES
            if (json_thispanel["nodes"].empty()) //TODO: Triangles
                    cout << "mull!" << endl;
            for (int index=0; index<4; index++)
                    thispanel->points[index] = *nodes + json_thispanel["nodes"][index].asInt();
            //is_wake
            thispanel->wake = json_thispanel["is_wake"].asBool();
            if (not thispanel->wake){
                    thispanel->position = j;
                    j++;
            }
    //set neighbours
            for (int index=0; index<4; index++){
        if (json_thispanel["neighbours"][index].isNull()){
            thispanel->neighbours[index] = NULL;
        } else {
            int num = json_thispanel["neighbours"][index].asInt();
            thispanel->neighbours[index] = *panels + num;
        }

            }
    }

    config->panel_number_not_wake = j;



    //CASES [v_inf1, v_inf2,...]
    //CONFIG: requests, farfield_factor, 



    return 1;



}

#endif