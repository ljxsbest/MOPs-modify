#include "swp_aluminum.h"

using namespace Sweep;
using namespace std;

//初始化列表
Aluminum::Aluminum():m_cap(0),r_cap(0){
	//no code.
}

//构造函数
Aluminum::Aluminum(void){
	cout << "constructor" << endl;
}

//析构函数
Aluminum::~Aluminum(void){
	cout << "destructor" << endl;
}

//设置初始质量
void set_init_mass(){

}

//设置初始环境温度（或许没必要）
void set_init_temp(){

}

//设置初始gas composition
void set_init_comp(){

}

//循环直到燃烧结束
for (unsigned int t = 0; t != endtime; ++t){

		calc_alum_mass(density,radius,t);
		calc_cap_mass(density,radius,t);

}