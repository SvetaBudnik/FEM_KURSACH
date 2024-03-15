#pragma once
#include <exception>
#include <vector>
using std::exception;

std::vector<double> sumVectors(const std::vector<double>& l, const std::vector<double>& r) {
   std::vector<double> result(l.size());
   for (int i = 0; i < l.size(); i++) {
      result[i] = l[i] + r[i];
   }
   return result;
}

std::vector<double> operator+ (const std::vector<double>& l, const std::vector<double>& r) {
   return sumVectors(l, r);
}

// Коэффициент диффузии (лямбда)
double diffusion(int omega, double r, double z) {
   double ans = 0;

   // выбираем функцию лямбды для соответствующей области
   switch (omega) {
   case 0: ans = r; break;

      // Если такой области омега нет, то
   default:
      throw std::exception("ТЫ ЗАБЫЛ/ЛА ПРО ЛЯМБДУ!!!");
   }
   return ans;
}


// Значение гаммы в области
//double gamma(int omega) {
//   double ans = 0;
//
//   // выбираем функцию гаммы для соответствующей области
//   switch (omega) {
//   case 0: ans = 1; break;
//
//      // Если такой области омега нет, то
//   default:
//      throw exception("ТЫ ЗАБЫЛ/ЛА ПРО gammu!!!");
//   }
//   return ans;
//}

// Значение кхи в области
double hi(int omega) {
   double ans = 0;

   // выбираем функцию хи для соответствующей области
   switch (omega) {
      case 0: ans = 0.001; break;

         // Если такой области омега нет, то
   default:
      throw exception("ТЫ ЗАБЫЛ/ЛА ПРО КХИ-КХИ!!!");
   }
   return ans;
}

// Значение сигма в области
double sigma(int omega) {
   double ans = 0;

   // выбираем функцию сыгма для соответствующей области
   switch (omega) {
      case 0: ans = 10; break;

         // Если такой области омега нет, то
   default:
      throw exception("ТЫ ЗАБЫЛ/ЛА ПРО СЫГМУ!!!");
   }
   return ans;
}


// Значение u_g в области (первые краевые)
double u_g(int omega, double r, double z, double t) {
   double ans = 0;

   // выбираем функцию уг для соответствующей области
   switch (omega) {
   case 0: ans = r + t; break;

      // Если такой области омега нет, то
   default:
      throw exception("ТЫ ЗАБЫЛ/ЛА ПРО U_G!!!!");
   }
   return ans;
}

double firstNachalnoe(int omega, double r, double z, double t) {
   double ans = 0;

   // выбираем функцию первого начального условия для соответствующей области
   switch (omega) {
   case 0: ans = u_g(omega, r, z, t); break;
   //case 0: ans = r * r + z * z; break;

      // Если такой области омега нет, то
   default:
      throw exception("ТЫ ЗАБЫЛ/ЛА ПРО первое начальное!!!!");
   }
   return ans;
}

//double secondNachalnoe(int omega, double r, double z, double t) {
//   double ans = 0;
//
//   // выбираем функцию второго начального условия (производной) для соответствующей области
//   switch (omega) {
//   //case 0: ans = r * r + z * z; break;
//
//      // Если такой области омега нет, то
//   default:
//      throw exception("ТЫ ЗАБЫЛ/ЛА ПРО второе начальное!!!!");
//   }
//   return ans;
//}

//
//// Значение бета в области (3 краевые) 
//double beta(int omega) {
//   double ans = 0;
//
//   // выбираем функцию бета для соответствующей области
//   switch (omega) {
//
//      // Если такой области омега нет, то
//   default:
//      throw exception("ТЫ ЗАБЫЛ/ЛА ПРО БЕТТУУУУУУ!!!");
//   }
//   return ans;
//}
//
//// Значение u_бета в области (3 краевые)
//double u_beta(int omega) {
//   double ans = 0;
//
//   // выбираем функцию бета для соответствующей области
//   switch (omega) {
//
//      // Если такой области омега нет, то
//   default:
//      throw exception("ТЫ ЗАБЫЛ/ЛА ПРО У_БЕТТУУУУУУ!!!");
//   }
//   return ans;
//}
//
//// Значение тета в области	(2 краевые)
//double tetta(int omega) {
//   double ans = 0;
//
//   // выбираем функцию тетта для соответствующей области
//   switch (omega) {
//
//      // Если такой области омега нет, то
//   default:
//      throw exception("ТЫ ЗАБЫЛ/ЛА ПРО ТАТЕТУУУ!!!");
//   }
//   return ans;
//}


// Значение функции f в точке
double f(int omega, double r, double z, double t) {
   double ans = 0;

   // выбираем функцию f для соответствующей области
   switch (omega) {
   case 0: ans = 8; break;

      // Если такой области омега нет, то
   default:
      throw exception("ТЫ ЗАБЫЛ/ЛА ПРО ФУНКЦИИ!!!");
   }
   return ans;
}
