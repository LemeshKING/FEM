#pragma once
#include "Header.h"
using namespace std;
class decision
{
public:
   void loadElements();
   void loadNodes();
   void loadBoundary();
   void loadTime();

private:
   int numberTime = 2;
   std::vector<Elemenets> elements;
   int amountElements;
   std::vector<std::vector<double>> nodes;
   int amountNodes;
   int amountTime;
   std::vector<std::vector<double>> A;
   std::vector<std::vector<double>> G;
   std::vector<std::vector<double>> Mass;
   std::vector<std::vector<double>> MassHi;
   std::vector<double> RightPart;
   vector<Elemenets> _mas;
   vector<vector<int>>_ltgm;
   std::vector<int> number_nodes;
   std::vector<int> ig;
   std::vector<int> jg;
   std::vector<double> gl;
   std::vector<double> di;
   std::vector<double> gu;
   std::vector<double> ggl;
   std::vector<double> ggu;
   std::vector<double> ddi;
   std::vector<double> GlobalRightPart;
   std::vector<std::vector<int>> boundary;
   std::vector<std::vector<int>> sBoundary;
   std::vector<std::vector<int>> ltgb;
   std::vector<double> localRightPartBoundary;
   std::vector < std::vector<double>> localMatrixBoundary;
   std::vector<double> x;
   std::vector<double> r0;
   std::vector<double> r;
   std::vector<double> z;
   std::vector<double> p;
   std::vector<double> help;
   std::vector<double> help2;
   std::vector<double> help3;
   std::vector<double> time;
   std::vector<double> Mdi; //j-1
   std::vector<double> Mgl; //j-1
   std::vector<double> Mgu; //j-1
   std::vector<double> Mdi1; //j-2
   std::vector<double> Mgl1; //j-2
   std::vector<double> Mgu1; //j-2
   std::vector<double> Mdi2; //j-3
   std::vector<double> Mgl2; //j-3
   std::vector<double> Mgu2; //j-3
   std::vector<std::vector<double>>qAll;
public:
   void print();
private:
   void getInitial();
   void getInitialForThree();
   void getInitialForFour();
   double function(int i, int reg, double time);
   double gamma(int i, int reg);
   double lambda(int reg);
   double Hi(int i, double time);
   void globalToLocal(int numberElement);
   void buildLocalGMatrix(int numberElement);
   void buildLocalMassMatrix(int numberElement);
   void buildLocalRightPart(int numberElement, double time);
   void buildLocalMatrix(int numberElement);
   void buildLocalMatrixHi();
   double theta(double x, double y);
   double u_g(double x, double y,double t);
   double ub(double x, double y);
public:
   void buildPortrait();
   void addLocalMatrix(int numberElement);
   void buildGlobalMatrix();
private:
   void addLocalMassMatrix(int numberElement, const std::vector<std::vector<double>> &Mq);
   void addLocalMassMatrix(int numberElement, const std::vector<std::vector<double>> &Mq, const std::vector<std::vector<double>> &Mq1);
   void addLocalMassMatrix(int numberElement, const std::vector<std::vector<double>> &Mq, const std::vector<std::vector<double>> &Mq1, const std::vector<std::vector<double>>& Mq2);
   void firstBoundaryCondition(double time);
   void secondBoundaryCondition();
   void thirdBoundaryCondition();
   void addLocalBoundary(int b);
public:
   void clear();
   void TwoLayerScheme();
   void ThreeLayerScheme();
   void solveLos();
   void profile2dense();
   void HyperbolicScheme();
   void FourLayerScheme();
private:
   void MqMult(const std::vector<double> &q, std::vector<double> &Mq); 
   void Mq1Mult(const std::vector<double> &q, std::vector<double> &Mq);
   void Mq2Mult(const std::vector<double>& q, std::vector<double>& Mq);
   void losInitialize();
   void lusqDecomposition();
   void solveL(std::vector<double> &x, const std::vector<double> &b);
   void solveU(std::vector<double> &x, const std::vector<double> &b);
   void multy(const std::vector<double> &x, std::vector<double> &result);
   double dotProduct(const std::vector<double> &x1, const std::vector<double> &x2);
   double sum(int i, int j);
};