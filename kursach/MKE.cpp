#include "MKE.h"
#include <fstream>
#include<set>
#include<algorithm>
#include<iostream>
#include<cmath>

using namespace std;

void decision::loadElements()
{
   ifstream input("elem.txt");
   input >> amountElements;
   elements.resize(amountElements);
   _mas.resize(amountElements);
   _ltgm.resize(amountElements);
   for (int i = 0; i < amountElements; i++)
   {
      elements[i].nodes.resize(6);
      _mas[i].nodes.resize(6);
      _ltgm[i].resize(6);
      for (int j = 0; j < 6; j++)
      {
         input >> elements[i].nodes[j];
         elements[i].nodes[j]--;
      }
      input >> elements[i].region;
      _mas[i] = elements[i];
      sort(_mas[i].nodes.begin(), _mas[i].nodes.begin() + 6);
      for (int j = 0; j < 6; j++) {
         int s;
         for (s = 0; s < 6; s++) {
            if (elements[i].nodes[s] == _mas[i].nodes[j])
            {
               _ltgm[i][j] = s;
               break;
            }
         }
      }
   }
   input.close();

  
}

void decision::loadNodes()
{
   ifstream input("nodes.txt");
   input >> amountNodes;
   nodes.resize(amountNodes);
   for (int i = 0; i < amountNodes; i++)
   {
      nodes[i].resize(2);
      input >> nodes[i][0] >> nodes[i][1];
   }
   input.close();
}

void decision::loadBoundary()
{
   ifstream input("boundary.txt");
  
   int n;
   input >> n;
   boundary.resize(n, vector<int>(5));
   sBoundary.resize(n, vector<int>(5));
   ltgb.resize(n, vector<int>(3));
   for (int i = 0; i < n; i++) 
   {
      for (int j = 0; j < 5; j++) 
      {
         input >> boundary[i][j];
         if (j > 1) --boundary[i][j];
      }
      sBoundary[i] = boundary[i];
      sort(sBoundary[i].begin() + 2, sBoundary[i].begin() + 5);
      for (int j = 2; j < 5; j++)
      {
         int s;
         for (s = 2; s < 5; s++)
            if (boundary[i][s] == sBoundary[i][j])
            {
               ltgb[i][j - 2] = s - 2;
               break;
            }
      }
   }
}

void decision::loadTime()
{
   ifstream input("time5.txt");
   input >> amountTime;
   time.resize(amountTime);
   qAll.resize(amountTime);
   for (int i = 0; i < amountTime; i++)
   {
      input >> time[i];
      qAll[i].resize(amountNodes);
   }
   
}

double decision::function(int i, int reg, double time)
{
   switch (reg)
   {
   case 1:
      return  nodes[i][0] * nodes[i][1];
         break;
   case 2: 
      return 1;
      break;
   case 3:
      return 1;
      break;
   default:
      break;
   }

}

double decision::gamma(int i, int reg)
{
   switch (reg)
   {
   case 1:
      return 1;
      break;
   case 2:
      return 0;
      break;
   case 3:
      return 0;
      break;
   default:
      break;
   }
}

double decision::lambda(int reg)
{
   switch (reg)
   {
   case 1:
      return 1;
      break;
   case 2:
      return 1;
      break;
   case 3:
      return 1;
      break;
   default:
      break;
   }
}

double decision::Hi(int i, double time)
{
   return 1;
}

void decision::globalToLocal(int numberElement)
{
   Elemenets sortElem = elements[numberElement];
   number_nodes.resize(6);
   sort(sortElem.nodes.begin(), sortElem.nodes.end());
   for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++)
         if (sortElem.nodes[i] == elements[numberElement].nodes[j])
            number_nodes[i] = elements[numberElement].nodes[i];
}

void decision::buildLocalGMatrix(int numberElement)
{
   int n1 = elements[numberElement].nodes[0];
   int n2 = elements[numberElement].nodes[1];
   int n3 = elements[numberElement].nodes[2];
   const double detD = (nodes[n2][0] - nodes[n1][0]) * (nodes[n3][1] - nodes[n1][1]) - (nodes[n3][0] - nodes[n1][0]) * (nodes[n2][1] - nodes[n1][1]);
   double A[3][3];
   double h1 = lambda(elements[numberElement].region);
   double h2 = lambda(elements[numberElement].region);
   double h3 = lambda(elements[numberElement].region);
   A[0][0] = nodes[n2][0] * nodes[n3][1] - nodes[n3][0] * nodes[n2][1];
   A[0][1] = nodes[n3][0] * nodes[n1][1] - nodes[n1][0] * nodes[n3][1];
   A[0][2] = nodes[n1][0] * nodes[n2][1] - nodes[n2][0] * nodes[n1][1];
   A[1][0] = nodes[n2][1] - nodes[n3][1]; A[2][0] = nodes[n3][0] - nodes[n2][0];
   A[1][1] = nodes[n3][1] - nodes[n1][1]; A[2][1] = nodes[n1][0] - nodes[n3][0];
   A[1][2] = nodes[n1][1] - nodes[n2][1]; A[2][2] = nodes[n2][0] - nodes[n1][0];
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
         A[i][j] /= detD;
   double a11 = (nodes[n2][1] - nodes[n3][1]) ;
   double a12 = (nodes[n3][0] - nodes[n2][0]);
   double a21 = (nodes[n3][1] - nodes[n1][1]) ;
   double a22 = (nodes[n1][0] - nodes[n3][0]) ;
   double a31 = (nodes[n1][1] - nodes[n2][1]) ;
   double a32 = (nodes[n2][0] - nodes[n1][0]) ;

   G.resize(6, vector<double>(6));
   G[0][4] = G[1][5] = G[2][3] = 0.0;
   G[0][0] = (A[1][0] * A[1][0] + A[2][0] * A[2][0]) / 2.0;
   G[0][1] = -(A[1][0] * A[1][1] + A[2][0] * A[2][1]) / 6.0;
   G[0][2] = -(A[1][0] * A[1][2] + A[2][0] * A[2][2]) / 6.0;
   G[0][3] = 2.0 * (A[1][0] * A[1][1] + A[2][0] * A[2][1]) / 3.0;
   G[0][5] = 2.0 * (A[1][0] * A[1][2] + A[2][0] * A[2][2]) / 3.0;
   G[1][1] = (A[1][1] * A[1][1] + A[2][1] * A[2][1]) / 2.0;
   G[1][2] = -(A[1][1] * A[1][2] + A[2][1] * A[2][2]) / 6.0;
   G[1][3] = 2.0 * (A[1][1] * A[1][0] + A[2][1] * A[2][0]) / 3.0;
   G[1][4] = 2.0 * (A[1][1] * A[1][2] + A[2][1] * A[2][2]) / 3.0;
   G[2][2] = (A[1][2] * A[1][2] + A[2][2] * A[2][2]) / 2.0;
   G[2][4] = 2.0 * (A[1][2] * A[1][1] + A[2][2] * A[2][1]) / 3.0;
   G[2][5] = 2.0 * (A[1][2] * A[1][0] + A[2][2] * A[2][0]) / 3.0;
   G[3][3] = 4.0 * (A[1][0] * A[1][0] + A[2][0] * A[2][0] + A[1][1] * A[1][1] + A[2][1] * A[2][1] + A[1][0] * A[1][1] + A[2][0] * A[2][1]) / 3.0;
   G[3][4] = 2.0 * (A[1][0] * A[1][1] + A[2][0] * A[2][1] + A[1][1] * A[1][2] + A[2][1] * A[2][2] + 2.0 * (A[1][0] * A[1][2] + A[2][0] * A[2][2]) + A[1][1] * A[1][1] + A[2][1] * A[2][1]) / 3.0;
   G[3][5] = 2.0 * (A[1][0] * A[1][0] + A[2][0] * A[2][0] + A[1][0] * A[1][2] + A[2][0] * A[2][2] + 2.0 * (A[1][1] * A[1][2] + A[2][1] * A[2][2]) + A[1][0] * A[1][1] + A[2][0] * A[2][1]) / 3.0;
   G[4][4] = 4.0 * (A[1][1] * A[1][1] + A[1][2] * A[1][2] + A[2][1] * A[2][1] + A[2][2] * A[2][2] + A[1][1] * A[1][2] + A[2][1] * A[2][2]) / 3.0;
   G[4][5] = 2.0 * (A[1][1] * A[1][2] + A[2][1] * A[2][2] + A[1][0] * A[1][2] + A[2][0] * A[2][2] + 2.0 * (A[1][1] * A[1][0] + A[2][1] * A[2][0]) + A[1][2] * A[1][2] + A[2][2] * A[2][2]) / 3.0;
   G[5][5] = 4.0 * (A[1][0] * A[1][0] + A[1][2] * A[1][2] + A[2][0] * A[2][0] + A[2][2] * A[2][2] + A[1][0] * A[1][2] + A[2][0] * A[2][2]) / 3.0;

   G[1][0] = G[0][1]; 
   G[2][0] = G[0][2];
   G[2][1] = G[1][2];
   G[3][0] = G[0][3]; 
   G[3][1] = G[1][3]; 
   G[3][2] = G[2][3];
   G[4][0] = G[0][4]; 
   G[4][1] = G[1][4]; 
   G[4][2] = G[2][4];
   G[4][3] = G[3][4];
   G[5][0] = G[0][5]; 
   G[5][1] = G[1][5];
   G[5][2] = G[2][5]; 
   G[5][3] = G[3][5]; 
   G[5][4] = G[4][5];
   for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++)
         G[i][j] *= h1 * abs(detD);
   


   for (int i = 0; i < 6; i++)
   
      cout << G[i][0] + G[i][1] + G[i][2] + G[i][3] + G[i][4] + G[i][5] << " ";
   cout << endl;
}

void decision::buildLocalMassMatrix(int numberElement)
{
   int reg = elements[numberElement].region;
   int n1 = elements[numberElement].nodes[0];
   int n2 = elements[numberElement].nodes[1];
   int n3 = elements[numberElement].nodes[2];
   const double detD = (nodes[n2][0] - nodes[n1][0]) * (nodes[n3][1] - nodes[n1][1]) - (nodes[n3][0] - nodes[n1][0]) * (nodes[n2][1] - nodes[n1][1]);
   double g1 = gamma(n1, reg), g2 = gamma(n2, reg), g3 = gamma(n3, reg);
   double tt = abs(detD);
   Mass.resize(6, vector<double>(6));
   //Mass[0][0] = abs(detD) * (5.0 * g1 + g2 + g3) / 420.0;
   //Mass[0][1] = abs(detD) * (g3 / 2520.0 - (g1 + g2) / 630.0);
   //Mass[0][2] = abs(detD) * (g2 / 2520.0 - (g1 + g3) / 630.0);
   //Mass[0][3] = abs(detD) * (g1 / 210.0 - g2 / 315 - g3 / 630.0);
   //Mass[0][4] = abs(detD) * (-g1 / 630.0 - (g2 + g3) / 210.0);
   //Mass[0][5] = abs(detD) * (g1 / 630.0 - (g2 + g3) / 210.0);
   //Mass[1][1] = abs(detD) * (5.0 * g2 + g1 + g3) / 420.0;
   //Mass[1][2] = abs(detD) * (g1 / 2520.0 - (g2 + g3) / 630.0);
   //Mass[1][3] = abs(detD) * (g2 / 210.0 - g1 / 315.0 - g3 / 630.0);
   //Mass[1][4] = abs(detD) * (g2 / 210 - g3 / 315.0 - g1 / 630);
   //Mass[1][5] = abs(detD) * (g2 / 630.0 - (g1 + g3) / 210.0);
   //Mass[2][2] = abs(detD) * (5 * g3 + g1 + g2) / 420.0;
   //Mass[2][3] = abs(detD) * (-g3 / 630.0 - (g1 + g2) / 210.0);
   //Mass[2][4] = abs(detD) * (g3 / 210.0 - g1 / 630.0 - g2 / 315.0);
   //Mass[2][5] = abs(detD) * (g3 / 210.0 - g1 / 315.0 - g2 / 630.0);
   //Mass[3][3] = abs(detD) * ((4 * g1 + 4 * g2) / 105.0 + 4 * g3 / 315.0);
   //Mass[3][4] = abs(detD) * ((4 * g1 + 4 * g3) / 315.0 + 2 * g2 / 105.0);
   //Mass[3][5] = abs(detD) * (2*g1/105.0 + (4*g2+4*g3)/315.0);
   //Mass[4][4] = abs(detD) * ((4 * g3 + 4 * g2) / 105 + 4 * g1 / 315.0);
   //Mass[4][5] = abs(detD) * (2 * g3 / 105.0 + (4 * g2 + 4 * g1) / 315.0);
   //Mass[5][5] = abs(detD) * ((4 * g1 + 4 * g3) / 105.0 + 4 * g2 / 315.0);
   //Mass[1][0] = Mass[0][1];
   //Mass[2][0] = Mass[0][2];
   //Mass[2][1] = Mass[1][2];
   //Mass[3][0] = Mass[0][3];
   //Mass[3][1] = Mass[1][3];
   //Mass[3][2] = Mass[3][2];
   //Mass[4][0] = Mass[0][4];
   //Mass[4][1] = Mass[1][4];
   //Mass[4][2] = Mass[2][4];
   //Mass[4][3] = Mass[3][4];
   //Mass[5][0] = Mass[0][5];
   //Mass[5][1] = Mass[1][5];
   //Mass[5][2] = Mass[2][5];
   //Mass[5][3] = Mass[3][5];
   //Mass[5][4] = Mass[4][5];
   Mass[0][0] = 1. / 60.;
   Mass[0][1] = -1. / 360.;
   Mass[0][2] = -1. / 360.;
   Mass[0][3] = 0;
   Mass[0][4] = -1. / 90.;
   Mass[0][5] = 0;
   Mass[1][0] = Mass[0][1];
   Mass[1][1] = 1. / 60.;
   Mass[1][2] = -1. / 360.;
   Mass[1][3] = 0;
   Mass[1][4] = 0;
   Mass[1][5] = -1. / 90;
   Mass[2][0] = Mass[0][2];
   Mass[2][1] = Mass[1][2];
   Mass[2][2] = 1. / 60;
   Mass[2][3] = -1. / 90.;
   Mass[2][4] = 0;
   Mass[2][5] = 0;
   Mass[3][0] = Mass[0][3];
   Mass[3][1] = Mass[1][3];
   Mass[3][2] = Mass[2][3];
   Mass[3][3] = 4. / 45.;
   Mass[3][4] = 2. / 45.;
   Mass[3][5] = 2. / 45.;
   Mass[4][0] = Mass[0][4];
   Mass[4][1] = Mass[1][4];
   Mass[4][2] = Mass[2][4];
   Mass[4][3] = Mass[3][4];
   Mass[4][4] = 4. / 45.;
   Mass[4][5] = 2. / 45.;
   Mass[5][0] = Mass[0][5];
   Mass[5][1] = Mass[1][5];
   Mass[5][2] = Mass[2][5];
   Mass[5][3] = Mass[3][5];
   Mass[5][4] = Mass[4][5];
   Mass[5][5] = 4. / 45.;
   MassHi.resize(6, vector<double>(6));
   for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++)
      {
         Mass[i][j] *= tt;
      }

}

void decision::buildLocalRightPart(int numberElement, double time)
{
   int reg = elements[numberElement].region;
   int n1 = elements[numberElement].nodes[0];
   int n2 = elements[numberElement].nodes[1];
   int n3 = elements[numberElement].nodes[2];
   const double detD = (nodes[n2][0] - nodes[n1][0]) * (nodes[n3][1] - nodes[n1][1]) - (nodes[n3][0] - nodes[n1][0]) * (nodes[n2][1] - nodes[n1][1]);
   
   double f1 = function(n1, reg, time), f2 = function(n2, reg, time), f3 = function(n3, reg, time),
      f4 = function(elements[numberElement].nodes[3], reg, time), f5 = function(elements[numberElement].nodes[4], reg, time), 
      f6 = function(elements[numberElement].nodes[5], reg, time);
   RightPart.resize(6);
   const double tt = abs(detD) / 360.;
  /* for (int i = 0; i < 6; i++)
      RightPart[i] = f1 * Mass[i][0] + f2 * Mass[i][1] + f3 * Mass[i][2] + f4 * Mass[i][3] + f5 * Mass[i][4] + f6 * Mass[i][5];*/
   RightPart[0] = tt * (6 * f1 - f2 - f3 - 4 * f5);
   RightPart[1] = tt * (6 * f2 - f1 - f3 - 4 * f6);
   RightPart[2] = tt * (6 * f3 - f1 - f2 - 4 * f4);
   RightPart[3] = tt * (32 * f4 + 16 * f5 + 16 * f6 - 4 * f3);
   RightPart[4] = tt * (32 * f5 + 16 * f4 + 16 * f6 - 4 * f1);
   RightPart[5] = tt * (32 * f6 + 16 * f4 + 16 * f5 - 4 * f2);
}

void decision::buildLocalMatrix(int numberElement)
{
   A.clear();
   A.resize(6, vector<double>(6));

   buildLocalGMatrix(numberElement);
   buildLocalMassMatrix(numberElement);
   for (int i = 0; i < 6; i++)
   {
      for (int j = 0; j < 6; j++)  
         A[i][j] += G[i][j] + Mass[i][j];
     
   }
   buildLocalRightPart(numberElement,0);
}

void decision::buildLocalMatrixHi()
{
   MassHi = Mass;
   for(int i = 0; i < 6; i++)
   {
      double hi = Hi(i, 0);
      for(int j = 0; j < 6; j++)
         MassHi[i][j] *= hi;
   }
   for (int i = 0; i < 6; i++)
   {
      double ksi = gamma(1,1);
      for (int j = 0; j < 6; j++)
         Mass[i][j] *= ksi;
   }
}

double decision::theta(double x, double y)
{
   return x+y;
}

double decision::u_g(double x, double y,double t)
{
   return x * y  * time[t];
}

void decision::print()
{
   vector<double> ideal;
   ofstream file("out.csv");
   file.precision(20);
   ideal.resize(amountNodes);
   for (int i = 0; i < amountTime; i++)
   {

         ideal[4] = u_g(nodes[4][0], nodes[4][1], i);
         file << nodes[4][0] << ";" << nodes[4][1] << ";";
         file << ideal[4] << ";" << qAll[i][4] << ";" << abs(qAll[i][4] - ideal[4]) << ";" << time[i];
         printf_s("%8.2E ", abs(qAll[i][4] - ideal[4]));
         cout << endl;

      cout << endl;
      file << endl;
   }
}

void decision::getInitial()
{
   for (int i = 0; i < amountNodes; i++)
      qAll[0][i] = u_g(nodes[i][0], nodes[i][1], 0);
}

void decision::getInitialForThree()
{
   for(int j = 0; j < 2; j++)
      for (int i = 0; i < amountNodes; i++)
         qAll[j][i] = u_g(nodes[i][0], nodes[i][1], j);
}

void decision::getInitialForFour()
{
   for(int j = 0; j < 3; j++)
      for(int i = 0; i < amountNodes; i++)
         qAll[j][i] = u_g(nodes[i][0],nodes[i][1],j);
}

double decision::ub(double x, double y)
{
   return x-y;
}

void decision::buildPortrait()
{

   vector<Elemenets> sortElem = elements;
   vector<set<int>> list(amountNodes);
   for (int i = 0; i < amountElements; i++)
   {
      sort(sortElem[i].nodes.begin(), sortElem[i].nodes.end());
      for (int j = 5; j >= 0; j--)
         for (int k = j - 1; k >= 0; k--)
            list[sortElem[i].nodes[j]].insert(sortElem[i].nodes[k]);
   }
   ig.resize(amountNodes + 1);
   ig[0] = 0;
   for (int i = 1; i < amountNodes + 1; i++)
      ig[i] = ig[i - 1] + list[i-1].size();
   Mdi.resize(amountNodes);
   Mgl.resize(ig[amountNodes]);
   Mgu.resize(ig[amountNodes]);
   jg.resize(ig[amountNodes]);
   di.resize(amountNodes);
   gl.resize(ig[amountNodes]);
   gu.resize(ig[amountNodes]);
   ggl.resize(ig[amountNodes]);
   ggu.resize(ig[amountNodes]);
   ddi.resize(amountNodes);
   GlobalRightPart.resize(amountNodes);
   int k = 0;
   for (int i = 0; i < amountNodes; i++)
      for (auto &elem : list[i])
      {
         jg[k] = elem;
         k++;
      }
}

void decision::addLocalMatrix(int numberElement)
{
   A.clear();
   A.resize(6, vector<double>(6));
   for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++)
         A[i][j] = G[i][j] + Mass[i][j] + MassHi[i][j];
   int k = 0;
   for (int i = 0; i < 6; i++)
   {
      di[elements[numberElement].nodes[i]] += A[i][i];
      GlobalRightPart[elements[numberElement].nodes[i]] += RightPart[i];
   }
   for (int i = 0; i < 6; i++)
      for (int j = 0; j < ig[elements[numberElement].nodes[i] + 1] - ig[elements[numberElement].nodes[i]]; j++)
      {
         bool flag = true;
         for (k = 0; k < 6 && flag; k++)
         {
            if (elements[numberElement].nodes[k] == jg[ig[elements[numberElement].nodes[i]] + j])
            {
               flag = false;
               break;
            }
         }
         if(!flag)
            gl[ig[elements[numberElement].nodes[i]] + j] += A[i][k];
      }

}

void decision::buildGlobalMatrix()
{
   for (int i = 0; i < amountElements; i++)
   {
      buildLocalMatrix(i);
      addLocalMatrix(i);
   }
   gu = gl;
   firstBoundaryCondition(0);
}

void decision::firstBoundaryCondition(double time)
{
   int i, n;
   double ug;
   for (i = 0, n = boundary.size(); i < n; ++i) 
      if (boundary[i][0] == 1) 
         break;
   const double BIG_NUM = 1;
   while (i < n && boundary[i][0] == 1) 
   {
      int region = boundary[i][1];
      for (int j = 2; j < 5; j++) 
      {
         ug = u_g(nodes[boundary[i][j]][0],nodes[boundary[i][j]][1],time);
         di[boundary[i][j]] = BIG_NUM;
         GlobalRightPart[boundary[i][j]] = ug;
         for (int k = 0; k < ig[boundary[i][j] + 1] - ig[boundary[i][j]]; k++)
            gl[ig[boundary[i][j]] + k] = 0;
         for (int k = 0; k < jg.size(); k++)
            if (boundary[i][j] == jg[k])
               gu[k] = 0;
      }
      ++i;
   }

}

void decision::secondBoundaryCondition()
{
   int i, n = boundary.size(), region;
   int n1 = 0, n2 = 0, n3 = 0;
   double tt = 0;
   localRightPartBoundary.resize(3);
   for (i = 0; i < n; ++i) 
      if (boundary[i][0] == 2)
         break;
   while (i < n && boundary[i][0] == 2)
   {
      region = boundary[i][1];
      n1 = boundary[i][2];
      n2 = boundary[i][3];
      n3 = boundary[i][4];
      
      double t1 = theta(nodes[n1][0],nodes[n1][1]), t2 = theta(nodes[n2][0], nodes[n2][1]), t3 = theta(nodes[n3][0], nodes[n3][1]);
       
      tt = sqrt((nodes[n1][0] - nodes[n3][0]) * (nodes[n1][0] - nodes[n3][0]) +
         (nodes[n1][1] - nodes[n3][1]) * (nodes[n1][1] - nodes[n3][1])) / 30.0;
      localRightPartBoundary[0] = tt * ((4) * t1 + (2) * t2 + (-1) * t3);
      localRightPartBoundary[1] = tt * ((2) * t1 + (16) * t2 + (2) * t3);
      localRightPartBoundary[2] = tt * ((-1) * t1 + (2) * t2 + (4) * t3);
      GlobalRightPart[n1] += localRightPartBoundary[0];
      GlobalRightPart[n2] += localRightPartBoundary[1];
      GlobalRightPart[n3] += localRightPartBoundary[2];
      ++i;
   }
}

void decision::thirdBoundaryCondition()
{
   int i, n = boundary.size(), region;
   int n1, n2, n3;
   double tt = 0;
   double beta = 1;
   localMatrixBoundary.resize(3, vector<double>(3));
   for (i = 0; i < n; ++i) 
      if (boundary[i][0] == 3)
         break;
   
   while (i < n && boundary[i][0] == 3)
   {
      
      region = boundary[i][1];
      n1 = boundary[i][2];
      n2 = boundary[i][3];
      n3 = boundary[i][4];
      double ub1 = ub(nodes[n1][0], nodes[n1][1]), ub2 = ub(nodes[n2][0], nodes[n2][1]), ub3 = ub(nodes[n3][0], nodes[n3][1]);
      tt = beta * sqrt((nodes[n1][0] - nodes[n3][0]) * (nodes[n1][0] - nodes[n3][0]) +
         (nodes[n1][1] - nodes[n3][1]) * (nodes[n1][1] - nodes[n3][1])) / 30.0;
      localRightPartBoundary[0] = tt * ((4) * ub1 + (2) * ub2 + (-1) * ub3);
      localRightPartBoundary[1] = tt * ((2) * ub1 + (16) * ub2 + (2) * ub3);
      localRightPartBoundary[2] = tt * ((-1) * ub1 + (2) * ub2 + (4) * ub3);
      GlobalRightPart[n1] += localRightPartBoundary[0];
      GlobalRightPart[n2] += localRightPartBoundary[1];
      GlobalRightPart[n3] += localRightPartBoundary[2];
      localMatrixBoundary[0][0] = tt * 4;
      localMatrixBoundary[0][1] = tt * 2;
      localMatrixBoundary[0][2] = tt * (-1);
      localMatrixBoundary[1][0] = tt * 2;
      localMatrixBoundary[1][1] = tt * 16;
      localMatrixBoundary[1][2] = tt * 2;
      localMatrixBoundary[2][0] = tt * (-1);
      localMatrixBoundary[2][1] = tt * 2;
      localMatrixBoundary[2][2] = tt * 4;
      addLocalBoundary(i);
      ++i;
   }
}

void decision::addLocalBoundary(int b)
{
   for (int i = 0; i < 3; ++i)
      di[boundary[b][i + 2]] += localMatrixBoundary[i][i];
   
      int j = 0;
      while (boundary[b][2] != jg[ig[boundary[b][3]] + j])
         j++;
      gl[ig[boundary[b][3]] + j] += localMatrixBoundary[1][0];
      int l = 0;
      while (boundary[b][2] != jg[ig[boundary[b][4]] + l])
         l++;
      gl[ig[boundary[b][4]] + l] += localMatrixBoundary[2][0];
      int m = 0;
      while (boundary[b][3] != jg[ig[boundary[b][4]] + m])
         m++;
      gl[ig[boundary[b][4]] + m] += localMatrixBoundary[2][1];
   
}

void decision::solveLos()
{
   lusqDecomposition();
   losInitialize();
   multy(x, help);
   for (int i = 0; i < amountNodes; ++i)
      r[i] = GlobalRightPart[i] - help[i];
   solveL(r0, r);
   r = r0;
   solveU(z, r);
   multy(z, help);
   solveL(p, help);
   double rr = 1;
   int k = 0;
   double eps = 1e-15;
   int maxIter = 10000000;
   eps *=eps;
   while (rr > eps && k < maxIter)
   {
      double pp = dotProduct(p, p);
      double alpha = dotProduct(p, r) / pp;
      for (int i = 0; i < amountNodes; ++i) {
         x[i] += alpha * z[i];
         r[i] -= alpha * p[i];
      }
      solveU(help, r);
      multy(help, help2);
      solveL(help3, help2);
      double beta = -dotProduct(p, help3) / pp;
      for (int i = 0; i < amountNodes; ++i) {
         z[i] = help[i] + beta * z[i];
         p[i] = help3[i] + beta * p[i];
      }
      if (k % 100 != 0)
         rr -= alpha * alpha * pp;
      else
         rr = dotProduct(r, r);
      ++k;
   }
   qAll[numberTime] = x;
   std::cout.precision(10);
   for (int i = 0; i < amountNodes; i++)
      std::cout << qAll[numberTime][i] << " ";
   numberTime++;
   std::cout << endl;
   std::cout << endl;
}

void decision::Mq1Mult(const std::vector<double> &q, std::vector<double> &Mq)
{
   int kk, i, j, kol_str, st;
   kk = 0;
   for (i = 0; i < amountNodes; i++)
      Mq[i] = Mdi1[i] * q[i];
   for (i = 1; i < amountNodes; i++) {
      kol_str = ig[i + 1] - ig[i];
      for (j = 0; j < kol_str; j++, kk++)
      {
         st = jg[kk];
         Mq[i] += Mgl1[kk] * q[st];
         Mq[st] += Mgu1[kk] * q[i];
      }
   }
}

void decision::Mq2Mult(const std::vector<double>& q, std::vector<double>& Mq)
{
   int kk, i, j, kol_str, st;
   kk = 0;
   for (i = 0; i < amountNodes; i++)
      Mq[i] = Mdi2[i] * q[i];
   for (i = 1; i < amountNodes; i++) {
      kol_str = ig[i + 1] - ig[i];
      for (j = 0; j < kol_str; j++, kk++)
      {
         st = jg[kk];
         Mq[i] += Mgl2[kk] * q[st];
         Mq[st] += Mgu2[kk] * q[i];
      }
   }
}

void decision::profile2dense()
{
   vector <vector <double>> A;
   A.resize(amountNodes); for (int i = 0; i < amountNodes; i++) A[i].resize(amountNodes);
   for (int i = 0; i < amountNodes; i++) for (int j = 0; j < amountNodes; j++) A[i][j] = 0.0;
   for (int i = 0; i < amountNodes; i++)
   {
      A[i][i] = di[i];
      for (int ii = ig[i]; ii < ig[i + 1]; ii++)
      {
         A[i][jg[ii]] = gl[ii];
         A[jg[ii]][i] = gu[ii];
      }
   }
   ofstream csv_out("GlobalMatrixPlot.csv");
   csv_out.precision(40);
   for (int i = 0; i < amountNodes; i++)
   {
      for (int j = 0; j < amountNodes; j++) csv_out << A[i][j] << ";";
      csv_out << endl;
   }
   csv_out << endl;
   csv_out.precision(40);
   for (int i = 0; i < amountNodes; i++) csv_out << GlobalRightPart[i] << ";";
   csv_out.close();
}

void decision::HyperbolicScheme()
{
   getInitialForThree();
   for (int i = 2; i < amountTime; i++)
   {
      clear();
      double deltaT = time[i] - time[i - 2];
      double deltaT1 = time[i - 1] - time[i - 2];
      double deltaT0 = time[i] - time[i - 1];
      double deltaTsq = deltaT0 * deltaT0;
      for (int j = 0; j < amountElements; j++)
      {
         buildLocalGMatrix(j);
         buildLocalMassMatrix(j);
         buildLocalMatrixHi();
         buildLocalRightPart(j, time[i]);
         vector<vector<double>>Mq1 = Mass;
         vector<vector<double>>Mq2 = Mass;
         for(int k = 0; k < 6; k++)
            for (int l = 0; l < 6; l++)
            {
               Mass[k][l] *= (deltaT + deltaT0) / (deltaT * deltaT0);
               Mq1[k][l] *= deltaT / (deltaT1 * deltaT0); // qj - 1
               Mq2[k][l] *= -deltaT0 / (deltaT * deltaT1) ; //qj - 2
               Mq1[k][l] += MassHi[k][l] * 2.0 / (deltaT1 * deltaT0); // qj - 1
               Mq2[k][l] += MassHi[k][l] * (- 2.0) / (deltaT1 * deltaT); //qj - 2
               MassHi[k][l] *= 2. / (deltaT * deltaT0);
            }
         addLocalMatrix(j);
         addLocalMassMatrix(j, Mq1, Mq2);
         gu = gl;
      }
      vector<double> Mq;
      Mq.resize(amountNodes);
      MqMult(qAll[i - 1], Mq);
      for (int k = 0; k < amountNodes; k++)
         GlobalRightPart[k] += Mq[k];
      Mq.clear();
      Mq.resize(amountNodes);
      Mq1Mult(qAll[i - 2], Mq);
      for (int k = 0; k < amountNodes; k++)
         GlobalRightPart[k] += Mq[k];
      firstBoundaryCondition(i);
      solveLos();
   }

}

void decision::FourLayerScheme()
{
   getInitialForFour();
   for (int i = 3; i < amountTime; i++)
   {
      clear();
      double deltaT0 = time[i] - time[i - 1];
      double deltaT1 = time[i] - time[i - 2];
      double deltaT2 = time[i] - time[i - 3];
      double deltaT3 = time[i - 1] - time[i - 2];
      double deltaT4 = time[i - 1] - time[i - 3];
      double deltaT5 = time[i - 2] - time[i - 3];
      for (int j = 0; j < amountElements; j++)
      {
         buildLocalGMatrix(j);
         buildLocalMassMatrix(j);
         buildLocalRightPart(j, time[i]);
         vector<vector<double>>Mq1 = Mass;
         vector<vector<double>>Mq2 = Mass;
         vector<vector<double>>Mq3 = Mass;
         for (int k = 0; k < 6; k++)
            for (int l = 0; l < 6; l++)
            {
               Mass[k][l] *= (deltaT0 * deltaT1 + deltaT0 * deltaT2 + deltaT1 * deltaT2) / (deltaT0 * deltaT1 * deltaT2);
               Mq1[k][l] *= (deltaT1 * deltaT2) / (deltaT0 * deltaT3 * deltaT4); // qj - 1
               Mq2[k][l] *= -(deltaT0 * deltaT2) / (deltaT1 * deltaT3 * deltaT5); //qj - 2
               Mq3[k][l] *= (deltaT0 * deltaT1) / (deltaT2 * deltaT4 * deltaT5); //qj - 3
            }
         addLocalMatrix(j);
         addLocalMassMatrix(j, Mq1, Mq2, Mq3);
         gu = gl;
      }
      vector<double> Mq;
      Mq.resize(amountNodes);
      MqMult(qAll[i - 1], Mq);
      for (int k = 0; k < amountNodes; k++)
         GlobalRightPart[k] += Mq[k];
      Mq.clear();
      Mq.resize(amountNodes);
      Mq1Mult(qAll[i - 2], Mq);
      for (int k = 0; k < amountNodes; k++)
         GlobalRightPart[k] += Mq[k];
      Mq.clear();
      Mq.resize(amountNodes);
      Mq2Mult(qAll[i - 3], Mq);
      for (int k = 0; k < amountNodes; k++)
         GlobalRightPart[k] += Mq[k];
      firstBoundaryCondition(i);
      solveLos();
   }
}

void decision::TwoLayerScheme()
{
   getInitial();
   //int t = 1;
   for (int t = 1; t < amountTime; t++)
   {
      clear();
      vector<double> Mq;
      Mq.resize(amountNodes);
      for (int i = 0; i < amountElements; i++)
      {


         int reg = elements[i].region;
         int n1 = elements[i].nodes[0];
         //int n2 = elements[i].nodes[1];
         //int n3 = elements[i].nodes[2];
         double g = gamma(n1, reg);
         buildLocalGMatrix(i);
         buildLocalMassMatrix(i);
         buildLocalRightPart(i, time[t]);
         double deltaT = time[t] - time[t - 1];
         for (int j = 0; j < 6; j++)
            for (int k = 0; k < 6; k++)
               Mass[j][k] *= 1.0 / deltaT;
         addLocalMassMatrix(i, Mass);
         addLocalMatrix(i);
         gu = gl;
      }
      MqMult(qAll[t - 1], Mq);
      for (int i = 0; i < amountNodes; i++)
         GlobalRightPart[i] += Mq[i];
      
      firstBoundaryCondition(t);
      profile2dense();
      solveLos();
   }
}

void decision::addLocalMassMatrix(int numberElement, const std::vector<std::vector<double>> &Mq, const std::vector<std::vector<double>> &Mq1)
{
   //int index[6];
   //vector<int> ll(6);
   //for (int k = 0; k < 6; k++)
   //{
   //   index[k] = _mas[numberElement].nodes[k];
   //   ll[k] = _ltgm[numberElement][k];
   //}
   //int i = 0;
   //for (i; i < 6; i++)
   //{
   //   Mdi[index[i]] += Mq[ll[i]][ll[i]];
   //   Mdi1[index[i]] += Mq1[ll[i]][ll[i]];
   //}
   //int j;
   //for (i = 0; i < 6; i++)
   //{
   //   int iBeg = ig[index[i]];
   //   for (j = 0; j < i; j++, iBeg++)
   //   {
   //      int iEnd = ig[index[i] + 1];
   //      while (jg[iBeg] != index[j])
   //      {
   //         int ind = (iBeg + iEnd) / 2;
   //         if (jg[ind] <= index[j])
   //            iBeg = ind;
   //         else
   //            iEnd = ind;
   //      }
   //      Mgl[iBeg] += Mq[ll[i]][ll[j]];
   //      Mgu[iBeg] += Mq[ll[j]][ll[i]];
   //      Mgl1[iBeg] += Mq1[ll[i]][ll[j]];
   //      Mgu1[iBeg] += Mq1[ll[j]][ll[i]];
   //   }
   //}
   int k = 0;
   for (int i = 0; i < 6; i++)
   {
      Mdi[elements[numberElement].nodes[i]] += Mq[i][i];
      Mdi1[elements[numberElement].nodes[i]] += Mq1[i][i];
   }
   for (int i = 0; i < 6; i++)
      for (int j = 0; j < ig[elements[numberElement].nodes[i] + 1] - ig[elements[numberElement].nodes[i]]; j++)
      {
         bool flag = true;
         for (k = 0; k < 6 && flag; k++)
         {
            if (elements[numberElement].nodes[k] == jg[ig[elements[numberElement].nodes[i]] + j])
            {
               flag = false;
               break;
            }
         }
         if (!flag)
         {
            Mgl[ig[elements[numberElement].nodes[i]] + j] += Mq[i][k];
            Mgu[ig[elements[numberElement].nodes[i]] + j] += Mq[k][i];
            Mgl1[ig[elements[numberElement].nodes[i]] + j] += Mq1[i][k];
            Mgu1[ig[elements[numberElement].nodes[i]] + j] += Mq1[k][i];
         }
      }

}

void decision::addLocalMassMatrix(int numberElement, const std::vector<std::vector<double>>& Mq, const std::vector<std::vector<double>>& Mq1, const std::vector<std::vector<double>>& Mq2)
{
   int k = 0;
   for (int i = 0; i < 6; i++)
   {
      Mdi[elements[numberElement].nodes[i]] += Mq[i][i];
      Mdi1[elements[numberElement].nodes[i]] += Mq1[i][i];
      Mdi2[elements[numberElement].nodes[i]] += Mq2[i][i];
   }
   for (int i = 0; i < 6; i++)
      for (int j = 0; j < ig[elements[numberElement].nodes[i] + 1] - ig[elements[numberElement].nodes[i]]; j++)
      {
         bool flag = true;
         for (k = 0; k < 6 && flag; k++)
         {
            if (elements[numberElement].nodes[k] == jg[ig[elements[numberElement].nodes[i]] + j])
            {
               flag = false;
               break;
            }
         }
         if (!flag)
         {
            Mgl[ig[elements[numberElement].nodes[i]] + j] += Mq[i][k];
            Mgu[ig[elements[numberElement].nodes[i]] + j] += Mq[k][i];
            Mgl1[ig[elements[numberElement].nodes[i]] + j] += Mq1[i][k];
            Mgu1[ig[elements[numberElement].nodes[i]] + j] += Mq1[k][i];
            Mgl2[ig[elements[numberElement].nodes[i]] + j] += Mq2[i][k];
            Mgu2[ig[elements[numberElement].nodes[i]] + j] += Mq2[k][i];
         }
      }
}


void decision::losInitialize()
{
   x.resize(amountNodes);
   r0.resize(amountNodes);
   r.resize(amountNodes);
   z.resize(amountNodes);
   p.resize(amountNodes);
   help.resize(amountNodes);
   help2.resize(amountNodes);
   help3.resize(amountNodes);
}

void decision::lusqDecomposition()
{
   ggl = gl;
   ggu = gu;
   ddi = di;
   int j;
   for (int i = 0; i < amountNodes; ++i) {
      for (j = ig[i]; j < ig[i + 1]; ++j)
      {
         gl[j] = (ggl[j] - sum(i, jg[j]));
         gu[j] = (ggu[j] - sum(jg[j], i)) / di[jg[j]];
      }
      di[i] = ddi[i] - sum(i, i);
   }
}

void decision::addLocalMassMatrix(int numberElement, const std::vector<std::vector<double>> &Mq)
{
   //int index[6];
   //vector<int> ll(6);
   //for (int k = 0; k < 6; k++)
   //{
   //   index[k] = _mas[numberElement].nodes[k];
   //   ll[k] = _ltgm[numberElement][k];
   //}
   //int i = 0;
   //for (i; i < 6; i++)
   //
   //   Mdi[index[i]] += Mq[ll[i]][ll[i]];

   //
   //int j;
   //for (i = 0; i < 6; i++)
   //{
   //   int iBeg = ig[index[i]];
   //   for (j = 0; j < i; j++, iBeg++)
   //   {
   //      int iEnd = ig[index[i] + 1];
   //      while (jg[iBeg] != index[j])
   //      {
   //         int ind = (iBeg + iEnd) / 2;
   //         if (jg[ind] <= index[j])
   //            iBeg = ind;
   //         else
   //            iEnd = ind;
   //      }
   //      Mgl[iBeg] += Mq[ll[i]][ll[j]];
   //      Mgu[iBeg] += Mq[ll[j]][ll[i]];

   //   }
   //}
   int k = 0;
   for (int i = 0; i < 6; i++)
   {
      Mdi[elements[numberElement].nodes[i]] += Mq[i][i];
   }
   for (int i = 0; i < 6; i++)
      for (int j = 0; j < ig[elements[numberElement].nodes[i] + 1] - ig[elements[numberElement].nodes[i]]; j++)
      {
         bool flag = true;
         for (k = 0; k < 6 && flag; k++)
         {
            if (elements[numberElement].nodes[k] == jg[ig[elements[numberElement].nodes[i]] + j])
            {
               flag = false;
               break;
            }
         }
         if (!flag)
         {
            Mgl[ig[elements[numberElement].nodes[i]] + j] += Mq[i][k];
            Mgu[ig[elements[numberElement].nodes[i]] + j] += Mq[k][i];
         }
      }
}



void decision::solveL(std::vector<double> &x, const std::vector<double> &b)
{
   int j;
   for (int i = 0; i < amountNodes; ++i) {
      double s = b[i];
      int jEnd = ig[i + 1];
      for (j = ig[i]; j < jEnd; ++j)
         s -= gl[j] * x[jg[j]];
      x[i] = s / di[i];
   }
}

void decision::solveU(std::vector<double> &x, const std::vector<double> &b)
{
   x = b;
   int j;
   for (int i = amountNodes - 1; i > -1; --i) {
      double xi = x[i];
      int jEnd = ig[i + 1];
      for (j = ig[i]; j < jEnd; ++j)
         x[jg[j]] -= gu[j] * xi;
   }
}

void decision::clear()
{
   di.clear();
   gl.clear();
   x.clear();
   r0.clear();
   r.clear();
   z.clear();
   p.clear();
   GlobalRightPart.clear();
   Mdi.clear();
   Mgl.clear();
   Mgu.clear();
   Mdi1.clear();
   Mgl1.clear();
   Mgu1.clear();
   Mdi2.clear();
   Mgl2.clear();
   Mgu2.clear();
   help.clear();
   help2.clear();
   help3.clear();
   di.resize(amountNodes);
   gl.resize(ig[amountNodes]);
   gu.resize(ig[amountNodes]);
   Mdi.resize(amountNodes);
   Mgl.resize(ig[amountNodes]);
   Mgu.resize(ig[amountNodes]);
   Mdi1.resize(amountNodes);
   Mgl1.resize(ig[amountNodes]);
   Mgu1.resize(ig[amountNodes]);
   Mdi2.resize(amountNodes);
   Mgl2.resize(ig[amountNodes]);
   Mgu2.resize(ig[amountNodes]);
   GlobalRightPart.resize(amountNodes);
}

void decision::multy(const std::vector<double> &x, std::vector<double> &result)
{
   int kk, i, j, kol_str, st;
   kk = 0;
   for (i = 0; i < amountNodes; i++)
      result[i] = ddi[i] * x[i];
   for (i = 1; i < amountNodes; i++) {
      kol_str = ig[i + 1] - ig[i];
      for (j = 0; j < kol_str; j++, kk++)
      {
         st = jg[kk];
         result[i] += ggl[kk] * x[st];
         result[st] += ggu[kk] * x[i];
      }
   }
}

double decision::dotProduct(const std::vector<double> &x1, const std::vector<double> &x2)
{
   double s = 0;
   for (int i = 0, n = x1.size(); i < n; ++i)
      s += x1[i] * x2[i];
   return s;
}

void decision::ThreeLayerScheme()
{
   getInitialForThree();
   for (int i = 2; i < amountTime; i++)
   {
      clear();
      double deltaT = time[i] - time[i - 2];
      double deltaT1 = time[i - 1] - time[i - 2];
      double deltaT0 = time[i] - time[i - 1];
      for (int j = 0; j < amountElements; j++)
      {
         buildLocalGMatrix(j);
         buildLocalMassMatrix(j);
         buildLocalRightPart(j, time[i]);
         vector<vector<double>> Mq1 = Mass;
         vector<vector<double>> Mq2 = Mass;
         for(int k = 0; k < 6; k++)
            for (int l = 0; l < 6; l++)
            {
               Mass[k][l] *= (deltaT + deltaT0) / (deltaT * deltaT0);
               Mq1[k][l] *= deltaT / (deltaT1 * deltaT0); //qj-1
               Mq2[k][l] *= -deltaT0 / (deltaT * deltaT1); //qj-2
            }
         addLocalMatrix(j);
         addLocalMassMatrix(j, Mq1, Mq2);
         gu = gl;

      }
      vector<double> Mq;
      Mq.resize(amountNodes);
      MqMult(qAll[i - 1], Mq);
      for (int k = 0; k < amountNodes; k++)
         GlobalRightPart[k] += Mq[k];
      Mq.clear();
      Mq.resize(amountNodes);
      Mq1Mult(qAll[i - 2], Mq);
      for (int k = 0; k < amountNodes; k++)
         GlobalRightPart[k] += Mq[k];
      firstBoundaryCondition(i);
      solveLos();
   }
}

double decision::sum(int i, int j)
{
   int k, l;
   bool find;
   double result = 0;
   if (i == j) {
      for (k = ig[i]; k < ig[i + 1]; ++k)
         result += gu[k] * gl[k];
   }
   else if (i > j) {
      // верхний треугольник
      for (k = ig[j]; k < ig[j + 1]; k++) {
         find = false;
         for (l = ig[i]; l < ig[i + 1]
            && !find; ++l) {
            if (jg[l] == jg[k]
               ) {
               result += gu[k] * gl[l];
               find = true;
            }
         }
      }
   }
   else {
      // нижний треугольник
      for (l = ig[i]; l < ig[i + 1]; l++)
      {
         find = false;
         for (k = ig[j]; k < ig[j + 1] && !find; ++k) 
         {
            if (jg[l] == jg[k])
            {
               result += gu[k] * gl[l];
               find = true;
            }
         }
      }
   }
   return result;
}

void decision::MqMult(const std::vector<double> &q, std::vector<double> &Mq)
{
   int kk, i, j, kol_str, st;
   kk = 0;
   for (i = 0; i < amountNodes; i++)
      Mq[i] = Mdi[i] * q[i];
   for (i = 1; i < amountNodes; i++) {
      kol_str = ig[i + 1] - ig[i];
      for (j = 0; j < kol_str; j++, kk++)
      {
         st = jg[kk];
         Mq[i] += Mgl[kk] * q[st];
         Mq[st] += Mgu[kk] * q[i];
      }
   }
}

