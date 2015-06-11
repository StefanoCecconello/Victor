/* 
 * File:   superImpositor.cc
 * Author: cecco
 * 
 * Created on June 8, 2015, 4:22 PM
 */

#include "SuperImpositor.h"


using namespace Victor;
using namespace Victor::Biopool;

SuperImpositor::SuperImpositor(Protein* firstProtein, Protein* secondProtein, string method="") {
    if (method == "Kabsch") {
        rotationAlgorith = new KabschMethod();
    }
    
    //set default rotation method
    if (method == "") {
        rotationAlgorith = new KabschMethod();
    }
    
    //Get the spacers
    set1 = firstProtein->getSpacer((unsigned int) 0);
    set2 = secondProtein->getSpacer((unsigned int) 0);
}

SuperImpositor::~SuperImpositor() {
}

double SuperImpositor::calculateRMSD(){
    
    
    
    
    Eigen::Affine3d* rotoTraslation=rotationAlgorith->rotate(set1,set2);
    rotoTraslation->linear();
    rotoTraslation->translation();
    return 0;
}

