#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

#include "IterSolvers.h" //Подключение ЛОС 
#include "gaussian_quadrature.h"
#include "functions.h"

using namespace std;

using Matrix = array<array<double, 6>, 6>;

//
//     1
//    / \
//   3   4
//  /     \
// 0 - 5 - 2
//
struct triangle {
   array<int, 6> inds;

   int omega = 0;
};

struct coordinates {
   double r = 0;
   double z = 0;
};

struct kraev_1 {
   int point = 0;
   int funcNum = 0;
};

//struct kraev_23 {
//   int point1 = 0;
//   int point2 = 0;
//   int funcNum = 0;
//};

vector <double> global_d;// Вектор правой части d
SparseMatrix global_matr; //Глобальная матрица в разреженном формате
vector<triangle> triangle_input; //Входные треугольники из файла triangle.txt
vector<coordinates> coord_input; //Входные данные координаты вершин point.txt (узлы)
vector<kraev_1> kraev1_input;
//vector<kraev_23> kraev2_input;
//vector<kraev_23> kraev3_input;

void input()
{
   ifstream fin;
   fin.open("triangles.txt");

   int size;// Количество треугольников
   fin >> size;//Считать размер и записать в size 

   triangle_input.resize(size);//Выделение места под все треугольники

   for (int i = 0; i < size; i++) { //Считывание треугольников
      for (int j = 0; j < triangle_input[i].inds.size(); j++) {
         fin >> triangle_input[i].inds[j];
      }
      fin >> triangle_input[i].omega;
   }
   fin.close();

   fin.open("point.txt");
   fin >> size;
   coord_input.resize(size);

   for (int i = 0; i < size; i++) {
      fin >> coord_input[i].r;
      fin >> coord_input[i].z;
   }

   fin.close();

   fin.open("kraev1.txt");
   fin >> size;
   kraev1_input.resize(size);

   for (int i = 0; i < size; i++) {
      fin >> kraev1_input[i].point;
      fin >> kraev1_input[i].funcNum;
   }
   fin.close();

   //fin.open("kraev2.txt");
   //fin >> size;
   //kraev2_input.resize(size);

   //for (int i = 0; i < size; i++) {
   //   fin >> kraev2_input[i].point1;
   //   fin >> kraev2_input[i].point2;
   //   fin >> kraev2_input[i].funcNum;
   //}
   //fin.close();

   //fin.open("kraev3.txt");
   //fin >> size;
   //kraev3_input.resize(size);

   //for (int i = 0; i < size; i++) {
   //   fin >> kraev3_input[i].point1;
   //   fin >> kraev3_input[i].point2;
   //   fin >> kraev3_input[i].funcNum;
   //}
   //fin.close();
}

// Генерация портрета глобальной матрицы
void generatePortrait() {
   global_matr.di.resize(coord_input.size());
   global_matr.ig.resize(coord_input.size() + 1);

   for (auto& triangle : triangle_input)
   {
      auto elems = triangle.inds;
      sort(begin(elems), end(elems));

      for (int i = 0; i < elems.size(); i++)
      {
         for (int k = 0; k < i; k++)
         {
            // Если элемент в верхнем прямоугольнике, то скипаем
            if (elems[k] > elems[i]) continue;

            bool isExist = false;
            // Пробегаем по всей строке для проверки, существует ли такой элемент
            for (auto it = global_matr.ig[elems[i]]; it < global_matr.ig[elems[i] + 1ll]; it++) {
               if (global_matr.jg[it] == elems[k]) {
                  isExist = true;
                  break;
               }
            }
            if (!isExist) {
               // Ищем, куда вставить элемент портрета
               auto it = global_matr.ig[elems[i]];
               while (it < global_matr.ig[elems[i] + 1ll] && global_matr.jg[it] < elems[k]) it++;

               // Для вставки нужно взять итератор массива от начала, так что...
               global_matr.jg.insert(global_matr.jg.begin() + it, elems[k]);

               // Добавляем всем элементам ig с позиции elems[i]+1 один элемент
               for (auto j = elems[i] + 1; j < global_matr.ig.size(); j++)
                  global_matr.ig[j]++;
            }
         }
      }
   }
   global_matr.ggl.resize(global_matr.jg.size());
   global_matr.ggu.resize(global_matr.jg.size());
}

double det_D(coordinates a, coordinates b, coordinates c)
{
   return ((b.r - a.r) * (c.z - a.z) - (c.r - a.r) * (b.z - a.z));
}

array<array<double, 3>, 4> getAlpha(const triangle& trig) {
   coordinates a = coord_input[trig.inds[0]];
   coordinates b = coord_input[trig.inds[1]];
   coordinates c = coord_input[trig.inds[2]];

   double det = det_D(a, b, c);

   array<array<double, 3>, 4> A;
   A[1][0] = (b.r * c.z - c.r * b.z) / det;
   A[1][1] = (b.z - c.z) / det;
   A[1][2] = (c.r - b.r) / det;
   A[2][0] = (c.r * a.z - a.r * c.z) / det;
   A[2][1] = (c.z - a.z) / det;
   A[2][2] = (a.r - c.r) / det;
   A[3][0] = (a.r * b.z - b.r * a.z) / det;
   A[3][1] = (a.z - b.z) / det;
   A[3][2] = (b.r - a.r) / det;

   return A;
}

array<double, 3> getL(const array<array<double, 3>, 4>& alphas, double r, double z) {
   array<double, 3> vec = { 1, r, z };
   array<double, 3> ans = {};
   for (int i = 1; i < 4; i++) {
      double sum = 0.0;
      for (int j = 0; j < 3; j++) {
         sum += alphas[i][j] * vec[j];
      }
      ans[i] = sum;
   }

   return ans;
}

double get_dPsi_dr(int ind, const triangle& trig, double r, double z) {
   auto alphas = getAlpha(trig);
   auto L = getL(alphas, r, z);

   if (ind == 0) {
      return 4 * alphas[0][1] * L[0] - alphas[0][1];
   }
   else if (ind == 1) {
      return 4 * alphas[1][1] * L[1] - alphas[1][1];
   }
   else if (ind == 2) {
      return 4 * alphas[2][1] * L[2] - alphas[2][1];
   }
   else if (ind == 3) {
      return 4 * (alphas[0][1] * L[1] + alphas[1][1] * L[0]);
   }
   else if (ind == 4) {
      return 4 * (alphas[1][1] * L[2] + alphas[2][1] * L[1]);
   }
   else {
      return 4 * (alphas[0][1] * L[2] + alphas[2][1] * L[0]);
   }
}

double get_dPsi_dz(int ind, const triangle& trig, double r, double z) {
   auto alphas = getAlpha(trig);
   auto L = getL(alphas, r, z);

   if (ind == 0) {
      return 4 * alphas[0][2] * L[0] - alphas[0][2];
   }
   else if (ind == 1) {
      return 4 * alphas[1][2] * L[1] - alphas[1][2];
   }
   else if (ind == 2) {
      return 4 * alphas[2][2] * L[2] - alphas[2][2];
   }
   else if (ind == 3) {
      return 4 * (alphas[0][2] * L[1] + alphas[1][2] * L[0]);
   }
   else if (ind == 4) {
      return 4 * (alphas[1][2] * L[2] + alphas[2][2] * L[1]);
   }
   else {
      return 4 * (alphas[0][2] * L[2] + alphas[2][2] * L[0]);
   }
}

double fact(int degree) {
   double res = 1.0;
   for (int i = 2; i <= degree; i++) {
      res *= i;
   }
   return res;
}

double integrateLLL(int c1, int c2, int c3, double det, double r1, double r2, double r3) {
   double res = 0.0;
   res += r1 * (fact(c1 + 1) * fact(c2) * fact(c3) / fact(c1 + c2 + c3 + 3));
   res += r2 * (fact(c1) * fact(c2 + 1) * fact(c3) / fact(c1 + c2 + c3 + 3));
   res += r3 * (fact(c1) * fact(c2) * fact(c3 + 1) / fact(c1 + c2 + c3 + 3));
   res *= abs(det);
   return res;
}

double getIntegral(int i, int j, const array<array<double, 3>, 4>& a, double det, triangle trig) {
   double r1 = coord_input[trig.inds[0]].r;
   double r2 = coord_input[trig.inds[1]].r;
   double r3 = coord_input[trig.inds[2]].r;

   auto lll = [&](int l1, int l2, int l3) {
      return integrateLLL(l1, l2, l3, det, r1, r2, r3);
      };

   double res = 0.0;

   if (i > j) {
      swap(i, j);
   }

   if (i == 0) {
      if (j == 0) {
         res += 16 * a[1][1] * a[1][1] * lll(2, 0, 0);
         res += 16 * a[1][2] * a[1][2] * lll(2, 0, 0);
         res -= 8 * a[1][1] * a[1][1] * lll(1, 0, 0);
         res -= 8 * a[1][2] * a[1][2] * lll(1, 0, 0);
         res += a[1][1] * a[1][1] * lll(0, 0, 0);
         res += a[1][2] * a[1][2] * lll(0, 0, 0);
         return res;
      }
      else if (j == 1) {
         res += 16 * a[1][1] * a[2][1] * lll(1, 1, 0);
         res += 16 * a[1][2] * a[2][2] * lll(1, 1, 0);
         res -= 4 * a[1][1] * a[2][1] * lll(1, 0, 0);
         res -= 4 * a[1][2] * a[2][2] * lll(1, 0, 0);
         res -= 4 * a[1][1] * a[2][1] * lll(0, 1, 0);
         res -= 4 * a[1][2] * a[2][2] * lll(0, 1, 0);
         res += a[1][1] * a[2][1] * lll(0, 0, 0);
         res += a[1][2] * a[2][2] * lll(0, 0, 0);
         return res;
      }
      else if (j == 2) {
         res += 16 * a[1][1] * a[3][1] * lll(1, 0, 1);
         res += 16 * a[1][2] * a[3][2] * lll(1, 0, 1);
         res -= 4 * a[1][1] * a[3][1] * lll(1, 0, 0);
         res -= 4 * a[1][2] * a[3][2] * lll(1, 0, 0);
         res -= 4 * a[1][1] * a[3][1] * lll(0, 0, 1);
         res -= 4 * a[1][2] * a[3][2] * lll(0, 0, 1);
         res += a[1][1] * a[3][1] * lll(0, 0, 0);
         res += a[1][2] * a[3][2] * lll(0, 0, 0);
         return res;
      }
      else if (j == 3) {
         res += 16 * a[1][1] * a[2][1] * lll(2, 0, 0);
         res += 16 * a[1][2] * a[2][2] * lll(2, 0, 0);
         res += 16 * a[1][1] * a[1][1] * lll(1, 1, 0);
         res += 16 * a[1][2] * a[1][2] * lll(1, 1, 0);
         res -= 4 * a[1][1] * a[2][1] * lll(1, 0, 0);
         res -= 4 * a[1][2] * a[2][2] * lll(1, 0, 0);
         res -= 4 * a[1][1] * a[1][1] * lll(0, 1, 0);
         res -= 4 * a[1][2] * a[1][2] * lll(0, 1, 0);
         return res;
      }
      else if (j == 4) {
         res += 16 * a[1][1] * a[3][1] * lll(1, 1, 0);
         res += 16 * a[1][2] * a[3][2] * lll(1, 1, 0);
         res += 16 * a[1][1] * a[2][1] * lll(1, 0, 1);
         res += 16 * a[1][2] * a[2][2] * lll(1, 0, 1);
         res -= 4 * a[1][1] * a[3][1] * lll(0, 1, 0);
         res -= 4 * a[1][2] * a[3][2] * lll(0, 1, 0);
         res -= 4 * a[1][1] * a[2][1] * lll(0, 0, 1);
         res -= 4 * a[1][2] * a[2][2] * lll(0, 0, 1);
         return res;
      }
      else if (j == 5) {
         res += 16 * a[1][1] * a[3][1] * lll(2, 0, 0);
         res += 16 * a[1][2] * a[3][2] * lll(2, 0, 0);
         res += 16 * a[1][1] * a[1][1] * lll(1, 0, 1);
         res += 16 * a[1][2] * a[1][2] * lll(1, 0, 1);
         res -= 4 * a[1][1] * a[3][1] * lll(1, 0, 0);
         res -= 4 * a[1][2] * a[3][2] * lll(1, 0, 0);
         res -= 4 * a[1][1] * a[1][1] * lll(0, 0, 1);
         res -= 4 * a[1][2] * a[1][2] * lll(0, 0, 1);
         return res;
      }
   }
   else if (i == 1) {
      if (j == 1) {
         res += 16 * a[2][1] * a[2][1] * lll(0, 2, 0);
         res += 16 * a[2][2] * a[2][2] * lll(0, 2, 0);
         res -= 8 * a[2][1] * a[2][1] * lll(0, 1, 0);
         res -= 8 * a[2][2] * a[2][2] * lll(0, 1, 0);
         res += a[2][1] * a[2][1] * lll(0, 0, 0);
         res += a[2][2] * a[2][2] * lll(0, 0, 0);
         return res;
      }
      else if (j == 2) {
         res += 16 * a[2][1] * a[3][1] * lll(0, 1, 1);
         res += 16 * a[2][2] * a[3][2] * lll(0, 1, 1);
         res -= 4 * a[2][1] * a[3][1] * lll(0, 1, 0);
         res -= 4 * a[2][2] * a[3][2] * lll(0, 1, 0);
         res -= 4 * a[2][1] * a[3][1] * lll(0, 0, 1);
         res -= 4 * a[2][2] * a[3][2] * lll(0, 0, 1);
         res += a[2][1] * a[3][1] * lll(0, 0, 0);
         res += a[2][2] * a[3][2] * lll(0, 0, 0);
         return res;
      }
      else if (j == 3) {
         res += 16 * a[2][1] * a[2][1] * lll(1, 1, 0);
         res += 16 * a[2][2] * a[2][2] * lll(1, 1, 0);
         res -= 4 * a[2][1] * a[2][1] * lll(1, 0, 0);
         res -= 4 * a[2][2] * a[2][2] * lll(1, 0, 0);
         res += 16 * a[1][1] * a[2][1] * lll(0, 2, 0);
         res += 16 * a[1][2] * a[2][2] * lll(0, 2, 0);
         res -= 4 * a[1][1] * a[2][1] * lll(0, 1, 0);
         res -= 4 * a[1][2] * a[2][2] * lll(0, 1, 0);
         return res;
      }
      else if (j == 4) {
         res += 16 * a[2][1] * a[3][1] * lll(0, 2, 0);
         res += 16 * a[2][2] * a[3][2] * lll(0, 2, 0);
         res += 16 * a[2][1] * a[2][1] * lll(0, 1, 1);
         res += 16 * a[2][2] * a[2][2] * lll(0, 1, 1);
         res -= 4 * a[2][1] * a[3][1] * lll(0, 1, 0);
         res -= 4 * a[2][2] * a[3][2] * lll(0, 1, 0);
         res -= 4 * a[2][1] * a[2][1] * lll(0, 0, 1);
         res -= 4 * a[2][2] * a[2][2] * lll(0, 0, 1);
         return res;
      }
      else if (j == 5) {
         res += 16 * a[2][1] * a[3][1] * lll(1, 1, 0);
         res += 16 * a[2][2] * a[3][2] * lll(1, 1, 0);
         res -= 4 * a[2][1] * a[3][1] * lll(1, 0, 0);
         res -= 4 * a[2][2] * a[3][2] * lll(1, 0, 0);
         res += 16 * a[1][1] * a[2][1] * lll(0, 1, 1);
         res += 16 * a[1][2] * a[2][2] * lll(0, 1, 1);
         res -= 4 * a[1][1] * a[2][1] * lll(0, 0, 1);
         res -= 4 * a[1][2] * a[2][2] * lll(0, 0, 1);
         return res;
      }
   }
   else if (i == 2) {
      if (j == 2) {
         res += 16 * a[3][1] * a[3][1] * lll(0, 0, 2);
         res += 16 * a[3][2] * a[3][2] * lll(0, 0, 2);
         res -= 8 * a[3][1] * a[3][1] * lll(0, 0, 1);
         res -= 8 * a[3][2] * a[3][2] * lll(0, 0, 1);
         res += a[3][1] * a[3][1] * lll(0, 0, 0);
         res += a[3][2] * a[3][2] * lll(0, 0, 0);
         return res;
      }
      else if (j == 3) {
         res += 16 * a[2][1] * a[3][1] * lll(1, 0, 1);
         res += 16 * a[2][2] * a[3][2] * lll(1, 0, 1);
         res -= 4 * a[2][1] * a[3][1] * lll(1, 0, 0);
         res -= 4 * a[2][2] * a[3][2] * lll(1, 0, 0);
         res += 16 * a[1][1] * a[3][1] * lll(0, 1, 1);
         res += 16 * a[1][2] * a[3][2] * lll(0, 1, 1);
         res -= 4 * a[1][1] * a[3][1] * lll(0, 1, 0);
         res -= 4 * a[1][2] * a[3][2] * lll(0, 1, 0);
         return res;
      }
      else if (j == 4) {
         res += 16 * a[3][1] * a[3][1] * lll(0, 1, 1);
         res += 16 * a[3][2] * a[3][2] * lll(0, 1, 1);
         res -= 4 * a[3][1] * a[3][1] * lll(0, 1, 0);
         res -= 4 * a[3][2] * a[3][2] * lll(0, 1, 0);
         res += 16 * a[2][1] * a[3][1] * lll(0, 0, 2);
         res += 16 * a[2][2] * a[3][2] * lll(0, 0, 2);
         res -= 4 * a[2][1] * a[3][1] * lll(0, 0, 1);
         res -= 4 * a[2][2] * a[3][2] * lll(0, 0, 1);
         return res;
      }
      else if (j == 5) {
         res += 16 * a[3][1] * a[3][1] * lll(1, 0, 1);
         res += 16 * a[3][2] * a[3][2] * lll(1, 0, 1);
         res -= 4 * a[3][1] * a[3][1] * lll(1, 0, 0);
         res -= 4 * a[3][2] * a[3][2] * lll(1, 0, 0);
         res += 16 * a[1][1] * a[3][1] * lll(0, 0, 2);
         res += 16 * a[1][2] * a[3][2] * lll(0, 0, 2);
         res -= 4 * a[1][1] * a[3][1] * lll(0, 0, 1);
         res -= 4 * a[1][2] * a[3][2] * lll(0, 0, 1);
         return res;
      }
   }
   else if (i == 3) {
      if (j == 3) {
         res += 16 * a[2][1] * a[2][1] * lll(2, 0, 0);
         res += 16 * a[2][2] * a[2][2] * lll(2, 0, 0);
         res += 32 * a[1][1] * a[2][1] * lll(1, 1, 0);
         res += 32 * a[1][2] * a[2][2] * lll(1, 1, 0);
         res += 16 * a[1][1] * a[1][1] * lll(0, 2, 0);
         res += 16 * a[1][2] * a[1][2] * lll(0, 2, 0);
         return res;
      }
      else if (j == 4) {
         res += 16 * a[2][1] * a[3][1] * lll(1, 1, 0);
         res += 16 * a[2][2] * a[3][2] * lll(1, 1, 0);
         res += 16 * a[2][1] * a[2][1] * lll(1, 0, 1);
         res += 16 * a[2][2] * a[2][2] * lll(1, 0, 1);
         res += 16 * a[1][1] * a[3][1] * lll(0, 2, 0);
         res += 16 * a[1][2] * a[3][2] * lll(0, 2, 0);
         res += 16 * a[1][1] * a[2][1] * lll(0, 1, 1);
         res += 16 * a[1][2] * a[2][2] * lll(0, 1, 1);
         return res;
      }
      else if (j == 5) {
         res += 16 * a[2][1] * a[3][1] * lll(2, 0, 0);
         res += 16 * a[2][2] * a[3][2] * lll(2, 0, 0);
         res += 16 * a[1][1] * a[3][1] * lll(1, 1, 0);
         res += 16 * a[1][2] * a[3][2] * lll(1, 1, 0);
         res += 16 * a[1][1] * a[2][1] * lll(1, 0, 1);
         res += 16 * a[1][2] * a[2][2] * lll(1, 0, 1);
         res += 16 * a[1][1] * a[1][1] * lll(0, 1, 1);
         res += 16 * a[1][2] * a[1][2] * lll(0, 1, 1);
         return res;
      }
   }
   else if (i == 4) {
      if (j == 4) {
         res += 16 * a[3][1] * a[3][1] * lll(0, 2, 0);
         res += 16 * a[3][2] * a[3][2] * lll(0, 2, 0);
         res += 32 * a[2][1] * a[3][1] * lll(0, 1, 1);
         res += 32 * a[2][2] * a[3][2] * lll(0, 1, 1);
         res += 16 * a[2][1] * a[2][1] * lll(0, 0, 2);
         res += 16 * a[2][2] * a[2][2] * lll(0, 0, 2);
         return res;
      }
      if (j == 5) {
         res += 16 * a[3][1] * a[3][1] * lll(1, 1, 0);
         res += 16 * a[3][2] * a[3][2] * lll(1, 1, 0);
         res += 16 * a[2][1] * a[3][1] * lll(1, 0, 1);
         res += 16 * a[2][2] * a[3][2] * lll(1, 0, 1);
         res += 16 * a[1][1] * a[3][1] * lll(0, 1, 1);
         res += 16 * a[1][2] * a[3][2] * lll(0, 1, 1);
         res += 16 * a[1][1] * a[2][1] * lll(0, 0, 2);
         res += 16 * a[1][2] * a[2][2] * lll(0, 0, 2);
         return res;
      }
   }
   else if (i == 5) {
      if (j == 5) {
         res += 16 * a[3][1] * a[3][1] * lll(2, 0, 0);
         res += 16 * a[3][2] * a[3][2] * lll(2, 0, 0);
         res += 32 * a[1][1] * a[3][1] * lll(1, 0, 1);
         res += 32 * a[1][2] * a[3][2] * lll(1, 0, 1);
         res += 16 * a[1][1] * a[1][1] * lll(0, 0, 2);
         res += 16 * a[1][2] * a[1][2] * lll(0, 0, 2);
         return res;
      }
   }

   throw runtime_error("Всякое случается, но такое....");
}

// Локальная матрица жёсткости
Matrix getLocalGFor(triangle loc_trig) {
   Matrix G;

   coordinates a = coord_input[loc_trig.inds[0]];
   coordinates b = coord_input[loc_trig.inds[1]];
   coordinates c = coord_input[loc_trig.inds[2]];

   int omega = loc_trig.omega;

   double diffusion_coef = diffusion(omega, a.r, a.z);

   double alpha[3][3];
   double det = abs(det_D(a, b, c));

   alpha[0][1] = (b.z - c.z) / det;
   alpha[0][2] = (c.r - b.r) / det;

   alpha[1][1] = (c.z - a.z) / det;
   alpha[1][2] = (a.r - c.r) / det;

   alpha[2][1] = (a.z - b.z) / det;
   alpha[2][2] = (b.r - a.r) / det;

   // TODO: реализация получения матрицы жёсткости численным интегрированием



#ifndef NDEBUG
   cout << "Проверка локальной матрицы жёсткости на корректность\n";
   for (int i = 0; i < 3; i++) {
      double tmp = 0.0;
      for (int j = 0; j < 3; j++) {
         tmp += G[i][j];
      }
      cout << "   " << i + 1 << ": сумма строки: " << tmp << endl;
   }
   cout << "\n";
#endif // !NDEBUG

   return G;
}

// Локальная матрица массы
Matrix getLocalM(triangle loc_trig, double gamma = 1.0) {
   Matrix m;
   coordinates a = coord_input[loc_trig.inds[0]];
   coordinates b = coord_input[loc_trig.inds[1]];
   coordinates c = coord_input[loc_trig.inds[2]];

   double g = hi(loc_trig.omega);

   double det = det_D(a, b, c);
   det = abs(det);

   // TODO: расчёт матрицы массы

   return m;
}

// Добавление локальной матрицы `local_mat` в глобальную `global_mat`
void addLocalToGlobal(triangle loc_trig, SparseMatrix& global_mat, Matrix& local_mat) {
   const auto& elems = loc_trig.inds;
   for (int i = 0; i < elems.size(); i++) {
      // добавляем все внедиагональные элементы на строке elems[i]
      for (int j = 0; j < i; j++) {
         int id;
         for (id = global_mat.ig[elems[i]]; id < global_mat.ig[elems[i] + 1] && global_mat.jg[id] != elems[j]; id++);

         global_mat.ggl[id] += local_mat[i][j];
         global_mat.ggu[id] += local_mat[j][i];
      }
      // добавляем диагональные элементы
      global_mat.di[elems[i]] += local_mat[i][i];
   }
}

// Локальный вектор правой части b
array<double, 3> getLocalbFor(triangle loc_trig, double t) {
   array<double, 3> vec_b;

   //coordinates a = coord_input[loc_trig.a];
   //coordinates b = coord_input[loc_trig.b];
   //coordinates c = coord_input[loc_trig.c];

   double f_loc[3] = {
      //f(loc_trig.omega, a.r, a.z, t),
      //f(loc_trig.omega, b.r, b.z, t),
      //f(loc_trig.omega, c.r, c.z, t),
   };

   //double det = abs(det_D(a, b, c));

   // TODO: сделать расчёт b через матрицу масс как b = Cf
   return vec_b;
}

// Первое краевое условие
void kr_1(double t) {
   for (int i = 0; i < kraev1_input.size(); i++) {
      int pos = kraev1_input[i].point;

      global_matr.di[pos] = 1;
      global_d[pos] = u_g(kraev1_input[i].funcNum, coord_input[pos].r, coord_input[pos].z, t);

      // зануляем строку в нижнем треугольнике
      for (int j = global_matr.ig[pos]; j < global_matr.ig[pos + 1]; j++)
      {
         global_matr.ggl[j] = 0;
      }
      // зануляем строку в верхнем треугольнике
      for (int i = pos + 1; i < global_matr.di.size(); i++)
      {
         for (int j = global_matr.ig[i]; j < global_matr.ig[i + 1]; j++)
         {
            if (global_matr.jg[j] == pos)
            {
               global_matr.ggu[j] = 0;
               break;
            }
         }
      }
   }
}

int main() {
   setlocale(LC_ALL, "ru-RU");

#ifndef NDEBUG
   cout << "********************************************" << endl;
   cout << "    программа запущена в режиме отладки     " << endl;
   cout << "********************************************" << endl;
   cout << "        *        " << endl;
   cout << "       ***       " << endl;
   cout << "      *   *      " << endl;
   cout << "      *   *      " << endl;
   cout << "      *   *      " << endl;
   cout << "      *   *      " << endl;
   cout << "      *   *      " << endl;
   cout << "      *   *      " << endl;
   cout << "   **** * ****   " << endl;
   cout << "   **   *   **   " << endl;
   cout << "  **    *    **  " << endl;
   cout << "     *     *     " << endl;
   cout << "    *       *    " << endl;
   cout << "   *         *   " << endl;
   cout << "Ввод данных начат..." << endl;
#endif // !NDEBUG

   // Делаем ввод данных и генерируем портрет
   input();
   generatePortrait();

   // Создаём временную сетку
   double t_begin = 0.0;
   double t_end = 3.0;
   double t_step = 1.0;
   vector<double> timeGrid((int)floor((t_end - t_begin) / t_step));
   for (int i = 0; i < timeGrid.size(); i++) {
      timeGrid[i] = t_begin + i * t_step;
   }

   // Создаём массив слоёв (кушек)
   vector<vector<double>> layers(3);

   // Находим первые два слоя (просто через начальные условия)
   layers[0].resize(coord_input.size());
   layers[1].resize(coord_input.size());
   for (int i = 0; i < coord_input.size(); i++) {
      auto& coords = coord_input[i];
      layers[0][i] = firstNachalnoe(0, coords.r, coords.z, timeGrid[0]);
      layers[1][i] = firstNachalnoe(0, coords.r, coords.z, timeGrid[1]);
   }

   // Находим глобальные матрицы
   SparseMatrix global_G = SparseMatrix::CopyShape(global_matr);
   SparseMatrix global_M_hi = SparseMatrix::CopyShape(global_matr);
   SparseMatrix global_M_sigma = SparseMatrix::CopyShape(global_matr);

   // цикл по всем треугольникам
   for (auto& triangle : triangle_input) {
      // Считаем локальную матрицу жёсткости
      Matrix local_g = getLocalGFor(triangle);
      addLocalToGlobal(triangle, global_G, local_g);

      // Считаем локальную матрицу массы хи
      Matrix local_m_hi = getLocalM(triangle, hi(triangle.omega));
      addLocalToGlobal(triangle, global_M_hi, local_m_hi);

      // Считаем локальную матрицу массы сигма
      Matrix local_m_sigma = getLocalM(triangle, sigma(triangle.omega));
      addLocalToGlobal(triangle, global_M_sigma, local_m_sigma);
   }

   //cout << "Глобальная жопа: \n" << global_G.toStringAsDense() << endl << endl;
   //cout << "Глобальная м хи: \n" << global_M_hi.toStringAsDense() << endl << endl;
   //cout << "Глобальная м сигма: \n" << global_M_sigma.toStringAsDense() << endl << endl;

   for (int j = 2; j < timeGrid.size(); j++) {
      cout << "Вычисления при t = " << timeGrid[j] << "\n\n";

      double t_j = timeGrid[j];
      double t_j_1 = timeGrid[j - 1];
      double t_j_2 = timeGrid[j - 2];

      double dt0 = t_j - t_j_1;
      double dt1 = t_j_1 - t_j_2;
      double dt2 = t_j - t_j_2;

      // Находим b_j
      vector<double> b_j(coord_input.size());
      for (auto& triangle : triangle_input) {
         array<double, 3> local_b_j = getLocalbFor(triangle, t_j);
         for (int i = 0; i < local_b_j.size(); i++) {
            b_j[triangle.inds[i]] += local_b_j[i];
         }
         //b_j[triangle.a] += local_b_j[0];
         //b_j[triangle.b] += local_b_j[1];
         //b_j[triangle.c] += local_b_j[2];
      }

      // Находим глобальную матрицу
      global_matr = global_M_hi * ((2.0) / (dt2 * dt0)) + global_M_sigma * ((dt2 + dt0) / (dt2 * dt0)) + global_G;
      //cout << "Глобальная матрица после сборки: \n" << global_matr.toStringAsDense() << "\n\n";

      // Находим вектор правой части

      //global_d = b_j + ((global_M_hi * (-2.0 / (dt1 * dt0))) * layers[0] + (global_M_hi * (2.0 / (dt1 * dt0))) * layers[1]) +
      //   ((global_M_sigma * (-dt2 / (dt1 * dt0))) * layers[0] + (global_M_sigma * (dt2 / (dt1 * dt0))) * layers[1]);
      global_d = b_j + ((global_M_hi * (-2.0 / (dt1 * dt2))) * layers[0] + (global_M_hi * (2.0 / (dt1 * dt0))) * layers[1]) +
         ((global_M_sigma * (-dt0 / (dt1 * dt2))) * layers[0] + (global_M_sigma * (dt2 / (dt1 * dt0))) * layers[1]);

      //cout << "Глобальный вектор после сборки: \n";
      //for (auto el : global_d) {
      //   cout << el << "\t";
      //}
      //cout << endl;

      // Накладываем первые краевые условия
      kr_1(t_j);
      //cout << "Глобальная матрица после наложения первых краевых:\n" << global_matr.toStringAsDense() << "\n\n";
      //cout << "Глобальный вектор после наложения краевых: \n";
      //for (auto el : global_d) {
      //   cout << el << "\t";
      //}
      //cout << endl;

      // Находим текущий слой
      layers[2] = vector<double>(b_j.size());
      double eps;
      IterSolvers::minEps = 1e-18;
      IterSolvers::maxIter = 500;
      IterSolvers::LOS::resetIter = 200;
      IterSolvers::LOS::Init_LuPrecond(b_j.size(), global_matr);
      IterSolvers::LOS::LuPrecond(global_matr, global_d, layers[2], eps);
      IterSolvers::Destruct();

      cout << "Результат работы программы: \n";
      for (int i = 0; i < layers[2].size(); i++) {
         cout << i << ": " << layers[2][i] << endl;
      }
      cout << endl << "Погрешность результата: \n";
      double nev = 0.0;
      for (int i = 0; i < layers[2].size(); i++) {
         nev += pow(layers[2][i] - u_g(0, coord_input[i].r, coord_input[i].z, t_j), 2);
         cout << i << ": " << layers[2][i] - u_g(0, coord_input[i].r, coord_input[i].z, t_j) << endl;
      }
      cout << endl << "Норма погрешности решения: " << sqrt(nev) << "\n\n\n";

      // Перемещаем слои
      layers[0] = layers[1];
      layers[1] = layers[2];
      // layers[2] занулится на следующей итерации
   }
}
