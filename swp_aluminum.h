#ifndef SWP_ALUMINUM_H_
#define SWP_ALUMINUM_H_

namespace Sweep{
	class Aluminum
	{
public:
	Aluminum(); //不含参构造函数声明
	~Aluminum(); //析构函数声明
	Aluminum(const Aluminum &copy); //拷贝构造函数
	//计算Al质量变化over time.
	void CalcAlMass();
	//计算cap质量变化over time.
	void CalcCapMass();
	//计算Al表面被cap包裹的ratio.
	double calcRatio();
		
	
	};

	
} //namespace Sweep

#endif