/**
 *  @file    Performances.cpp
 *  @brief   Define performances tests.
 *  @author  Francois Roy
 *  @date    12/01/2019
 */
#include "Performances.hpp"

namespace bench{
namespace performances{

Eigen::VectorXd loop(){
    Eigen::VectorXd d = Eigen::VectorXd::Zero(n);
    for (int i = 0; i < n; i++) {
        d[i] = u[i+1] - u[i];       
    }
    return d;
}
Eigen::VectorXd vectorized_loop(){
	// Eigen::VectorXd u = Eigen::VectorXd::Constant(n + 1, 2.4);
	Eigen::VectorXd d = Eigen::VectorXd::Zero(n);
	d = u.segment(1, u.size()-1) - u.segment(0, u.size()-1);
    return d;
}
std::vector<int> left_v(int nx, int ny){
    std::vector<int> temp; 
    for (int j=0; j<ny + 1; ++j){
        temp.push_back(j * (nx + 1));
    }
    return temp;
}
std::vector<int> bottom_v(int nx, int ny){
    std::vector<int> temp; 
    for (int i=0; i<nx + 1; ++i){
        temp.push_back(i);
    }
    return temp;
}
int bottom_left_corner_v(int nx, int ny){
	std::vector<int> left = left_v(nx, ny);
	std::vector<int> bottom = bottom_v(nx, ny);
    // need to sort first, maybe a better to do that...
	std::sort(left.begin(), left.end());
    std::sort(bottom.begin(), bottom.end());
    std::vector<int> s3;
    std::set_intersection(left.begin(),left.end(),bottom.begin(),bottom.end(), 
    	std::back_inserter(s3));
    return s3[0];
}
int bottom_left_corner_v_1(int nx, int ny){
	// only sort the smaller vector. Then do a single pass over the 
    // bigger vector and test a presence of its items in a smaller vector 
    // by using a binary search.
    std::vector<int> small, large, s3;
    if(nx <= ny){
    	small = left_v(nx, ny);
	    large = bottom_v(nx, ny);
    } else{
    	large = left_v(nx, ny);
	    small = bottom_v(nx, ny);
    }
	std::sort(small.begin(), small.end());
	for(int i=0; i<large.size(); i++){
		if (std::binary_search (small.begin(), small.end(), large[i])){
			s3.push_back(large[i]);
		}
	}
    return s3[0];
}
std::vector<int> union_v(int nx, int ny){
	std::vector<int> left = left_v(nx, ny);
	std::vector<int> bottom = bottom_v(nx, ny);
    // need to sort first, maybe a better to do that...
	std::sort(left.begin(), left.end());
    std::sort(bottom.begin(), bottom.end());
    std::vector<int> s3;
    std::set_union(left.begin(),left.end(),bottom.begin(),bottom.end(), 
    	std::back_inserter(s3));
    return s3;
}
std::set<int> left_s(int nx, int ny){
    std::set<int> temp; 
    for (int j=0; j<ny + 1; ++j){
        temp.insert(j * (nx + 1));
    }
    return temp;
}
std::set<int> bottom_s(int nx, int ny){
    std::set<int> temp; 
    for (int i=0; i<nx + 1; ++i){
        temp.insert(i);
    }
    return temp;
}
int bottom_left_corner_s(int nx, int ny){
	std::set<int> left = left_s(nx, ny);
	std::set<int> bottom = bottom_s(nx, ny);
    std::vector<int> s3;
    std::set_intersection(left.begin(),left.end(),bottom.begin(),bottom.end(), 
    	std::back_inserter(s3));
    return s3[0];
}
std::vector<int> union_s(int nx, int ny){
	std::set<int> left = left_s(nx, ny);
	std::set<int> bottom = bottom_s(nx, ny);
    std::vector<int> s3;
    std::set_union(left.begin(),left.end(),bottom.begin(),bottom.end(), 
    	std::back_inserter(s3));
    return s3;
}
Eigen::VectorXd uleft_s(int nx, int ny){
	int n = (nx+1)*(ny+1);
	Eigen::VectorXd u = Eigen::VectorXd::LinSpaced(n,0,n-1);
	std::set<int> s = left_s(nx, ny);
    std::set<int>::iterator setIt = s.begin();
    for( auto it = s.begin(); it!=s.end(); ++it){
        int ans = *it;
        u[ans] = 999.9;
     }
    return u;
}
Eigen::VectorXd uleft_v(int nx, int ny){
	int n = (nx+1)*(ny+1);
	Eigen::VectorXd u = Eigen::VectorXd::LinSpaced(n,0,n-1);
	std::vector<int> v = left_v(nx, ny);
    for( int i = 0; i<v.size(); ++i){
        u[i] = 999.9;
     }
    return u;
}

}  // namespace performances
}  // namespace bench
