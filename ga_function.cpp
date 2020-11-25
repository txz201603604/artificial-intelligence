#include "pch.h"
//=============================���ò���====================================
void ga::SetParameters()
{
	N_genetic = 50;                //��Ⱥ��ģ��̫С������̬������Ⱥ��ģ̫������������һ��0-100
	M_pgentic = 0.25;            //������ʣ�����Ⱥ�������йأ�һ��0.0001-0.2
	C_pgentic = 0.5;             //������ʣ�����̫�����״�ʧ���Ÿ��壬̫С�����Ǹ���Ч������Ⱥ��һ��0.4-0.99.
	E_gentic = 400;                 //����������̫С���㷨������������̫������ʱ�����Դ�˷ѣ�һ��100-500.
	precision = 0.001;
	L_variable = int(log((x_high[0] - x_low[0]) / precision + 1) / log(2));//����������ַ�������(������)
	N_variable = 2;
}
//***************************
//������ת��Ϊʵ��
//******************************
vector<double> ga::Real_trans(vector<int>x_binary)
{
	vector<int>x_decimal(N_variable);
	vector<double>x(N_variable);
	for (int j = 0; j < N_variable; j++)
	{
		for (int k = j * L_variable, l_gen = 0; k < (j + 1) * L_variable; k++, l_gen++)
		{
			x_decimal[j] = x_decimal[j] + x_binary[k] * pow(2, l_gen);
		}
		x[j] = x_low[0] + double(x_decimal[j]) / (pow(2, L_variable) - 1) * (x_high[0] - x_low[0]);
	}
	return x;
}
//*******************
//��ʼ������ֵ
//*******************
void ga::initialize()
{
	extern RandomNumber r;       //����ȫ�������
	x_i.resize(N_genetic, vector<double>(N_variable));
	x_best.resize(N_variable);
	fitness.resize(N_genetic);
	x_binary.resize(N_genetic, vector<int>(N_variable * L_variable)); //�Ż�����������
	for (int i = 0; i < N_genetic; i++)
	{
		//================================�������===============================================
		for (int j = 0; j < N_variable * L_variable; j++)
		{
			x_binary[i][j] = r.integer(0, 1);
			cout << x_binary[i][j];
		}
		cout << endl;
		x_i[i] = Real_trans(x_binary[i]);
		fitness[i] = 1 / function(x_i[i]);
	}
	x_best = x_i[0];                      //��ʼ�����Ÿ���
	best_fitness = fitness[0];
	for (int i = 1; i < N_genetic; i++)
	{
		if (best_fitness < fitness[i])
		{
			best_fitness = fitness[i];
			x_best = x_i[i];
		}
	}
}
void ga::Optimization_iteration()
{
	clock_t startTime, endTime; //�������ʼ����ʱ��ͽ���ʱ��
	startTime = clock();        //��ʱ��ʼ
	ofstream out("�Ŵ��㷨�Ż����.txt");
	for (int i = 0; i < E_gentic; i++)
	{
		select_operator();  //ѡ�񸸴����¸���
		crossover_operator();  //�н���Ȩ�����и������н���
		mutate_operator();     //�������
		for (int j = 0; j < N_genetic; j++)
		{
			fitness[j] = 1 / function(x_i[j]);
			if (best_fitness < fitness[j])
			{
				best_fitness = fitness[j];
				x_best = x_i[j];
			}
		}   //��Ⱥѡ��
		out << i << fixed << setw(12) << setprecision(5) << 1 / best_fitness << endl;

	}
	out << "���ű���:" << endl;
	for (int i = 0; i < N_variable; i++)
	{
		out << "x" << i << "=" << fixed << setw(12) << setprecision(5) << x_best[i] << endl;//������ű���
	}
	out << "����ֵ=" << fixed << setw(12) << setprecision(5) << 1 / best_fitness << endl;
	endTime = clock();//��ʱ����
	out << "run time:" << (double)(endTime - startTime) / CLOCKS_PER_SEC << "s" << endl;
	out.close();
}

//*******************
//ѡ�����Ӻ���
//*******************
void ga::select_operator()
{

	double totalfit = 0, p_i;
	sumfitness.resize(N_genetic);
	P_i.resize(N_genetic);
	extern RandomNumber r;       //�����
	vector<vector<int>>new_x_binary(N_genetic, vector<int>(N_variable * L_variable));//����ѡ��������¸���
	for (int i = 0; i < N_genetic; i++)
	{
		sumfitness[i] = totalfit + fitness[i];
		totalfit = totalfit + fitness[i];
	}
	//����������
	for (int i = 0; i < N_genetic; i++)
	{
		P_i[i] = sumfitness[i] / totalfit;
	}
	//ѡ�񸸴�
	for (int i = 0; i < N_genetic; i++)
	{
		p_i = r.decimal(0, 1.0);
		//�������̷�ѡ�����
		if (p_i <= P_i[0])
			new_x_binary[i] = x_binary[0];
		else
		{
			for (int j = 0; j < N_genetic - 1; j++)
			{
				if (p_i > P_i[j] && p_i <= P_i[j + 1])
					new_x_binary[i] = x_binary[j + 1];
			}
		}
	}
	//���¸���
	x_binary = new_x_binary;
}
//*******************
//�������Ӻ��������㽻��
//*******************
void ga::crossover_operator()
{
	int  cpoint1, cpoint2, t;                //�����cpoint1, cpoint2����,tΪ�滻ֵ
	double p_c;                              //��������������
	extern RandomNumber r;       //�����
	for (int i = 0; i < N_genetic; i = i + 2)
	{//�������������������
		cpoint1 = r.integer(0, N_variable * L_variable - 1);
		cpoint2 = r.integer(0, N_variable * L_variable - 1);
		if (cpoint2 < cpoint1)
		{
			t = cpoint2; cpoint2 = cpoint1; cpoint1 = t;
		}
		p_c = r.decimal(0, 1.0);
		//�������
		if (p_c < C_pgentic)
		{
			for (int j = cpoint1; j <= cpoint2; j++)
			{
				t = x_binary[i][j]; x_binary[i][j] = x_binary[i + 1][j]; x_binary[i + 1][j] = t;
			}
		}
	}
}
//*******************
//�������Ӻ�
//*******************
void ga::mutate_operator()
{
	int  mpoint;//�����mpoint����
	double p_m;
	extern RandomNumber r;       //����ȫ�������
	for (int i = 0; i < N_genetic; ++i)
	{//������������
		mpoint = r.integer(0, N_variable * L_variable - 1);
		p_m = r.decimal(0, 1.0);
		if (p_m < M_pgentic)
			if (x_binary[i][mpoint] == 0)
			{
				x_binary[i][mpoint] = 1;
			}
			else {
				x_binary[i][mpoint] = 0;
			}
	}
	//������Ⱦɫת��Ϊʵ��
	for (int i = 0; i < N_genetic; i++)
	{
		x_i[i] = Real_trans(x_binary[i]);
	}
}