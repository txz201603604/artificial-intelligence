#include "pch.h"
int main()
{
	ga GA;          //����ȫ���Ŵ��㷨��غ���
	//=========================�����㷨����====================
	GA.SetParameters();
	//========================��ʼ������ֵ=====================
	GA.initialize();
	//=====================�Ż�������������=======================
	GA.Optimization_iteration();
}