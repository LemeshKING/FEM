#include<iostream>
#include "MKE.h"

using namespace std;


int main()
{
   decision decision;
   decision.loadElements();
   decision.loadNodes();
   decision.loadBoundary();
   decision.loadTime();
   decision.buildPortrait();
   //decision.buildGlobalMatrix();
   //decision.solveLos();
   //decision.TwoLayerScheme();
   //decision.ThreeLayerScheme();
   decision.HyperbolicScheme();
   //decision.FourLayerScheme();
   decision.print();
   decision.profile2dense();
   return 0;
}