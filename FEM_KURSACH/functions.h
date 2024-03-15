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

// ����������� �������� (������)
double diffusion(int omega, double r, double z) {
   double ans = 0;

   // �������� ������� ������ ��� ��������������� �������
   switch (omega) {
   case 0: ans = r; break;

      // ���� ����� ������� ����� ���, ��
   default:
      throw std::exception("�� �����/�� ��� ������!!!");
   }
   return ans;
}


// �������� ����� � �������
//double gamma(int omega) {
//   double ans = 0;
//
//   // �������� ������� ����� ��� ��������������� �������
//   switch (omega) {
//   case 0: ans = 1; break;
//
//      // ���� ����� ������� ����� ���, ��
//   default:
//      throw exception("�� �����/�� ��� gammu!!!");
//   }
//   return ans;
//}

// �������� ��� � �������
double hi(int omega) {
   double ans = 0;

   // �������� ������� �� ��� ��������������� �������
   switch (omega) {
      case 0: ans = 0.001; break;

         // ���� ����� ������� ����� ���, ��
   default:
      throw exception("�� �����/�� ��� ���-���!!!");
   }
   return ans;
}

// �������� ����� � �������
double sigma(int omega) {
   double ans = 0;

   // �������� ������� ����� ��� ��������������� �������
   switch (omega) {
      case 0: ans = 10; break;

         // ���� ����� ������� ����� ���, ��
   default:
      throw exception("�� �����/�� ��� �����!!!");
   }
   return ans;
}


// �������� u_g � ������� (������ �������)
double u_g(int omega, double r, double z, double t) {
   double ans = 0;

   // �������� ������� �� ��� ��������������� �������
   switch (omega) {
   case 0: ans = r + t; break;

      // ���� ����� ������� ����� ���, ��
   default:
      throw exception("�� �����/�� ��� U_G!!!!");
   }
   return ans;
}

double firstNachalnoe(int omega, double r, double z, double t) {
   double ans = 0;

   // �������� ������� ������� ���������� ������� ��� ��������������� �������
   switch (omega) {
   case 0: ans = u_g(omega, r, z, t); break;
   //case 0: ans = r * r + z * z; break;

      // ���� ����� ������� ����� ���, ��
   default:
      throw exception("�� �����/�� ��� ������ ���������!!!!");
   }
   return ans;
}

//double secondNachalnoe(int omega, double r, double z, double t) {
//   double ans = 0;
//
//   // �������� ������� ������� ���������� ������� (�����������) ��� ��������������� �������
//   switch (omega) {
//   //case 0: ans = r * r + z * z; break;
//
//      // ���� ����� ������� ����� ���, ��
//   default:
//      throw exception("�� �����/�� ��� ������ ���������!!!!");
//   }
//   return ans;
//}

//
//// �������� ���� � ������� (3 �������) 
//double beta(int omega) {
//   double ans = 0;
//
//   // �������� ������� ���� ��� ��������������� �������
//   switch (omega) {
//
//      // ���� ����� ������� ����� ���, ��
//   default:
//      throw exception("�� �����/�� ��� ����������!!!");
//   }
//   return ans;
//}
//
//// �������� u_���� � ������� (3 �������)
//double u_beta(int omega) {
//   double ans = 0;
//
//   // �������� ������� ���� ��� ��������������� �������
//   switch (omega) {
//
//      // ���� ����� ������� ����� ���, ��
//   default:
//      throw exception("�� �����/�� ��� �_����������!!!");
//   }
//   return ans;
//}
//
//// �������� ���� � �������	(2 �������)
//double tetta(int omega) {
//   double ans = 0;
//
//   // �������� ������� ����� ��� ��������������� �������
//   switch (omega) {
//
//      // ���� ����� ������� ����� ���, ��
//   default:
//      throw exception("�� �����/�� ��� ��������!!!");
//   }
//   return ans;
//}


// �������� ������� f � �����
double f(int omega, double r, double z, double t) {
   double ans = 0;

   // �������� ������� f ��� ��������������� �������
   switch (omega) {
   case 0: ans = 8; break;

      // ���� ����� ������� ����� ���, ��
   default:
      throw exception("�� �����/�� ��� �������!!!");
   }
   return ans;
}
