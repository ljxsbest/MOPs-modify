#ifndef SWP_ALUMINUM_H_
#define SWP_ALUMINUM_H_

namespace Sweep{
	class Aluminum
	{
	public:
		Aluminum(); //不含参构造函数声明
		~Aluminum(); //析构函数声明
		Aluminum(const Aluminum &copy); //拷贝构造函数
	
		Aluminum(void){
		cout << "object is being created" << endl;
	}


	};
//成员函数定义，构造函数。
	Aluminum::Aluminum(void){
		cout << "object is being created" << endl;
	}
	//计算质量变化。
	void Aluminum::calc_mass(){

	}
	//计算初始cap的radius，及其所占比例belta.
	double Aluminum::calc_radius(a,b){

		return Rc, Belta
	}
	
}

#endif