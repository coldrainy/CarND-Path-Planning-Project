#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

//add function
vector<string> successor_states(string state,int lane) {
    /*
    Provides the possible next states given the current state for the FSM
    discussed in the course, with the exception that lane changes happen
    instantaneously, so LCL and LCR can only transition back to KL.
    */
    vector<string> states;
    states.push_back("KL");
    if(state.compare("KL") == 0) {
        states.push_back("PLCL");
        states.push_back("PLCR");
    } else if (state.compare("PLCL") == 0) {
        if (lane != 0) {
            states.push_back("PLCL");
            states.push_back("LCL");
        }
    } else if (state.compare("PLCR") == 0) {
        if (lane != 2) {
            states.push_back("PLCR");
            states.push_back("LCR");
        }
    } else if(state.compare("LCL") == 0){
        if (lane != 0)
            states.push_back("LCL");
    } else if(state.compare("LCR") == 0){
        if (lane != 2)
            states.push_back("LCR");
    }
    //If state is "LCL" or "LCR", then just return "KL"
    return states;
}
struct point{
    double x;
    double y;
};

map<string, int> lane_direction = {{"KL", 0} , {"PLCL", -1}, {"LCL", -1},
                                    {"LCR", 1}, {"PLCR", 1}};

bool get_vehicle_ahead(vector<vector<double>>sensor_fusion,int lane,vector<double> car_info, vector<double>& car_ahead_info,int previous_size){
    int min_s = 7000;
    bool found_vehicle = false;
    for(auto car:sensor_fusion){
        double check_d = car[6];
        if(check_d>(2+4*lane-2) && check_d<(2+4*lane+2) ){
            double check_s = car[5];
            double check_vx = car[3];
            double check_vy = car[4];
            double check_v = sqrt(check_vx*check_vx+check_vy*check_vy);
            check_s += check_v*0.02*previous_size;
            if(check_s > car_info[6] && check_s < min_s){
                min_s = check_s;
                car_ahead_info = car;
                found_vehicle = true;
            }
        }
    }
    return found_vehicle;
}
bool get_vehicle_behind(vector<vector<double>>sensor_fusion,int lane,vector<double> car_info, vector<double>& car_behind_info,int previous_size){
    int max_s = -1;
    bool found_vehicle = false;
    for(auto car:sensor_fusion){
        double check_d = car[6];
        if (check_d>(2+4*lane-2) && check_d<(2+4*lane+2)) {
            double check_s = car[5];
            double check_vx = car[3];
            double check_vy = car[4];
            double check_v = sqrt(check_vx*check_vx+check_vy*check_vy);
            check_s += check_v*0.02*previous_size;
            if(check_s < car_info[6] && check_s > max_s)
                max_s = check_s;
                car_behind_info = car;
                found_vehicle = true;
        }
    }
    return found_vehicle;
}
void get_kinematics(vector<double> car_info,vector<vector<double>> sensor_fusion,
        int lane,double& ref_val,int previous_size) {
    double max_velocity_accel_limit = ref_val + 5*0.02*2.24;
    double car_end_d = car_info[5];
    vector<double> car_ahead_info;
    vector<double> car_behind_info;
    if(get_vehicle_ahead(sensor_fusion, lane, car_info, car_ahead_info,previous_size)){
        double car_s = car_info[2];
        double car_end_s = car_info[6];
        if(previous_size>0)
            car_s = car_end_s;
        double car_ahead_s = car_ahead_info[5];
        double car_ahead_vx = car_ahead_info[3];
        double car_ahead_vy = car_ahead_info[4];
        double car_ahead_v = sqrt(car_ahead_vx*car_ahead_vx+car_ahead_vy*car_ahead_vy);
        double check_s = car_ahead_s+car_ahead_v*0.02*previous_size;
        double max_velocity_in_front = ref_val;
        if(check_s > car_s && check_s - car_s < 30){
            std::cout<<"speed down"<<std::endl;
            if(max_velocity_in_front > car_ahead_v )
                max_velocity_in_front -= 0.224;
        }
        else if(ref_val < 49){
            std::cout<<"speed up"<<std::endl;
            max_velocity_in_front += 0.224;
        }
        ref_val = min(min(max_velocity_in_front, max_velocity_accel_limit), (double)49);

        //if(get_vehicle_behind(sensor_fusion, lane, car_info, car_behind_info,previous_size)){
        if(0){

            double car_behind_s = car_behind_info[5];
            double car_behind_vx = car_behind_info[3];
            double car_behind_vy = car_behind_info[4];
            double car_behind_v = sqrt(car_behind_vx*car_behind_vx+car_behind_vy*car_behind_vy);
            double check_behind_s = car_behind_s+car_behind_v*0.02*previous_size;
            ref_val = car_ahead_v;
        }
    }else{
        std::cout<<"no vehicle ahead"<<std::endl;
        ref_val = min(max_velocity_accel_limit,(double)49);
    }
}

vector<point> getTrajectory(int lane,vector<vector<double>> sensor_fusion,
                            vector<double> car_info,vector<double> previous_path_x,
                            vector<double> previous_path_y,
                            double end_path_s,double end_path_d,
                            vector<double> map_waypoints_s,
                            vector<double> map_waypoints_x,
                            vector<double> map_waypoints_y,
                            double ref_val){
    vector<point> trajectory;

    vector<double> ptsx;
    vector<double> ptsy;

    double car_x = car_info[0];
    double car_y = car_info[1];
    double car_s = car_info[2];
    double car_d = car_info[3];
    double car_yaw = car_info[4];
    double car_speed = car_info[5];

    double ref_x = car_x;
    double ref_y = car_y;
    double ref_yaw = deg2rad(car_yaw);

    int previous_size = previous_path_x.size();
    //if(previous_size > 0){
        //car_s = end_path_s;
    //}
    //bool too_close = false;
    //for(auto obj:sensor_fusion){
        //double check_d = obj[6];
        //if( check_d>(2+4*lane-2) && check_d<(2+4*lane+2)){
            //double check_s = obj[5];
            //double check_vx = obj[3];
            //double check_vy = obj[4];
            //double check_v = sqrt(check_vx*check_vx + check_vy*check_vy);

            //check_s += check_v*0.02*previous_size;

            //if(check_s > car_s && check_s - car_s < 30){
                ////ref_val = 29.5;
                //too_close = true;
                ////if(lane>0)
                    ////lane = 0;
            //}
        //}
    //}
    //if(too_close){
        //ref_val -= 0.224;
    //}else if(ref_val < 49){
        //ref_val +=0.224;
    //}
    if(previous_size < 2){
        double pre_car_x = car_x - cos(ref_yaw);
        double pre_car_y = car_y - sin(ref_yaw);

        ptsx.push_back(pre_car_x);
        ptsx.push_back(car_x);

        ptsy.push_back(pre_car_y);
        ptsy.push_back(car_y);

    }else{
        ref_x = previous_path_x[previous_size-1];
        ref_y = previous_path_y[previous_size-1];

        double pre_ref_x = previous_path_x[previous_size-2];
        double pre_ref_y = previous_path_y[previous_size-2];

        ref_yaw = atan2(ref_y - pre_ref_y,ref_x - pre_ref_x);

        ptsx.push_back(pre_ref_x);
        ptsx.push_back(ref_x);

        ptsy.push_back(pre_ref_y);
        ptsy.push_back(ref_y);
    }
    if(previous_size > 0)
        car_s = car_info[6];
    vector<double> next_wp0 = getXY(car_s+30,2+4*lane,map_waypoints_s,map_waypoints_x,map_waypoints_y);
    vector<double> next_wp1 = getXY(car_s+60,2+4*lane,map_waypoints_s,map_waypoints_x,map_waypoints_y);
    vector<double> next_wp2 = getXY(car_s+90,2+4*lane,map_waypoints_s,map_waypoints_x,map_waypoints_y);

    ptsx.push_back(next_wp0[0]);
    ptsx.push_back(next_wp1[0]);
    ptsx.push_back(next_wp2[0]);

    ptsy.push_back(next_wp0[1]);
    ptsy.push_back(next_wp1[1]);
    ptsy.push_back(next_wp2[1]);

    for(int i=0;i<ptsx.size();i++){
        double shift_x = ptsx[i] - ref_x;
        double shift_y = ptsy[i] - ref_y;

        ptsx[i] = (shift_x*cos(0-ref_yaw) - shift_y*sin(0-ref_yaw));
        ptsy[i] = (shift_x*sin(0-ref_yaw) + shift_y*cos(0-ref_yaw));
    }
    tk::spline s;

    s.set_points(ptsx,ptsy);
    for(int i=0;i<previous_size;i++){
        trajectory.push_back(point{previous_path_x[i],previous_path_y[i]});
        //next_x_vals.push_back(previous_path_x[i]);
        //next_y_vals.push_back(previous_path_y[i]);
    }

    double target_x = 30;
    double target_y = s(target_x);
    double dist = sqrt(target_x*target_x+target_y*target_y);

    double x_add_on = 0;
    for(int i=0;i<50-previous_size;i++){

        int N = dist/(0.02*ref_val/2.24);
        double x_point = x_add_on + target_x/N;
        double y_point = s(x_point);

        x_add_on = x_point;

        double x_ref = x_point;
        double y_ref = y_point;

        x_point = (x_ref*cos(ref_yaw) - y_ref*sin(ref_yaw));
        y_point = (x_ref*sin(ref_yaw) + y_ref*cos(ref_yaw));

        x_point += ref_x;
        y_point += ref_y;

        trajectory.push_back(point{x_point,y_point});
        //next_x_vals.push_back(x_point);
        //next_y_vals.push_back(y_point);
    }
    return trajectory;

}
vector<point> keep_lane_trajectory(int lane,string state,vector<vector<double>> sensor_fusion,
                                   vector<double> car_info,vector<double> previous_path_x,
                                   vector<double> previous_path_y,
                                   double end_path_s,double end_path_d,
                                   vector<double> map_waypoints_s,
                                   vector<double> map_waypoints_x,
                                   vector<double> map_waypoints_y,
                                   double& ref_val,vector<float>& intend_lane_info,
                                   vector<float>& final_lane_info){
    std::cout<<"before kinematic ref_val:"<<ref_val<<std::endl;
    //double car_end_d = car_info[7];
    //int car_end_lane = car_end_d/4;
    get_kinematics(car_info,sensor_fusion,lane,ref_val,previous_path_x.size());
    intend_lane_info.push_back(ref_val);
    intend_lane_info.push_back(lane);
    final_lane_info.push_back(ref_val);
    final_lane_info.push_back(lane);

    std::cout<<"after kinematic ref_val:"<<ref_val<<std::endl;
    return getTrajectory(lane,sensor_fusion,car_info,previous_path_x,
                  previous_path_y,end_path_s,end_path_d,map_waypoints_s,
                  map_waypoints_x,map_waypoints_y,ref_val);
}
vector<point> lane_change_trajectory(int lane,string state,vector<vector<double>> sensor_fusion,
                                   vector<double> car_info,vector<double> previous_path_x,
                                   vector<double> previous_path_y,
                                   double end_path_s,double end_path_d,
                                   vector<double> map_waypoints_s,
                                   vector<double> map_waypoints_x,
                                   vector<double> map_waypoints_y,
                                   double& ref_val,vector<float>& intend_lane_info,vector<float>& final_lane_info){
    vector<point> trajectory;
    //double car_end_d = car_info[7];
    //int car_end_lane = car_end_d/4;

    int next_lane = lane + lane_direction[state];
    double car_end_s = car_info[6];
    double car_s = car_info[2];
    for(auto vehicle:sensor_fusion){
        double d = vehicle[6];
        double s = vehicle[5];
        double vx = vehicle[3];
        double vy = vehicle[4];
        double v = sqrt(vx*vx+vy*vy);
        double check_s = s+v*0.02*previous_path_x.size();


        //if(int(d/4) == next_lane && s == car_s)
        //if(d>(2+4*next_lane-2) && d<(2+4*lane+2) && check_s > car_end_s-20 && check_s < car_end_s+20){
        if(d>(2+4*next_lane-2) && d<(2+4*next_lane+2)){
            std::cout<<"d: "<<d<<"s: "<<s<<"car_s: "<<car_s<<std::endl;
            if(s > car_s-10 && s < car_s+5){
                intend_lane_info.push_back(0);
                intend_lane_info.push_back(-1);
                final_lane_info.push_back(0);
                final_lane_info.push_back(-1);

                std::cout<<"there is obstacle in the lane to change"<<std::endl;
                return trajectory;
            }
        }
    }
    std::cout<<"before kinematic ref_val:"<<ref_val<<std::endl;
    get_kinematics(car_info,sensor_fusion,next_lane,ref_val,previous_path_x.size());
    std::cout<<"after kinematic ref_val:"<<ref_val<<std::endl;
    intend_lane_info.push_back(ref_val);
    intend_lane_info.push_back(next_lane);
    final_lane_info.push_back(ref_val);
    final_lane_info.push_back(next_lane);
    return getTrajectory(next_lane,sensor_fusion,car_info,previous_path_x,
                  previous_path_y,end_path_s,end_path_d,map_waypoints_s,
                  map_waypoints_x,map_waypoints_y,ref_val);
}
vector<point> prep_lane_change_trajectory(int lane,string state,vector<vector<double>> sensor_fusion,
                                   vector<double> car_info,vector<double> previous_path_x,
                                   vector<double> previous_path_y,
                                   double end_path_s,double end_path_d,
                                   vector<double> map_waypoints_s,
                                   vector<double> map_waypoints_x,
                                   vector<double> map_waypoints_y,
                                   double& ref_val,vector<float>& intend_lane_info,vector<float>& final_lane_info){

    std::cout<<"before kinematic ref_val:"<<ref_val<<std::endl;
    //double car_end_d = car_info[7];
    //int car_end_lane = car_end_d/4;
    int next_lane = lane + lane_direction[state];

    double car_end_lane_val = ref_val;
    get_kinematics(car_info,sensor_fusion,lane,car_end_lane_val,previous_path_x.size());
    std::cout<<"after kinematic this_lane_val:"<<car_end_lane_val<<std::endl;
    vector<double> car_behind_info;
    //if(get_vehicle_behind(sensor_fusion,lane,car_info,car_behind_info,previous_size)){
    if(0){
        ref_val = car_end_lane_val;
    }else{
        double next_lane_val = ref_val;
        get_kinematics(car_info,sensor_fusion,next_lane,next_lane_val,previous_path_x.size());
        std::cout<<"after kinematic next_lane_val:"<<next_lane_val<<std::endl;
        if(car_end_lane_val > next_lane_val){
            ref_val = next_lane_val;
        }else{
            ref_val = car_end_lane_val;
        }
        intend_lane_info.push_back(next_lane_val);
        intend_lane_info.push_back(next_lane);
        final_lane_info.push_back(car_end_lane_val);
        final_lane_info.push_back(lane);
    }
    return getTrajectory(lane,sensor_fusion,car_info,previous_path_x,
                  previous_path_y,end_path_s,end_path_d,map_waypoints_s,
                  map_waypoints_x,map_waypoints_y,ref_val);
}
vector<point> generate_trajectory(int lane,string state,
                                  vector<vector<double>> sensor_fusion,
                                  vector<double> car_info,
                                  vector<double> previous_path_x,
                                  vector<double> previous_path_y,
                                  double end_path_s,
                                  double end_path_d,
                                  vector<double> map_waypoints_s,
                                  vector<double> map_waypoints_x,
                                  vector<double> map_waypoints_y,
                                  double& ref_val,vector<float> &intend_lane_info,
                                  vector<float>& final_lane_info){

    vector<point> trajectory;
    //int intended_lane = lane + lane_direction[state];
    if(state == "KL"){
        std::cout<<"------------------line 460: "<<state<<"------------------------"<<std::endl;
        trajectory = keep_lane_trajectory(lane,state,sensor_fusion,car_info,
                                  previous_path_x,previous_path_y,end_path_s,end_path_d,
                                  map_waypoints_s,map_waypoints_x,map_waypoints_y,ref_val,
                                  intend_lane_info,final_lane_info);
    }else if(state == "LCL" || state == "LCR"){
        std::cout<<"------------------line 466: "<<state<<" ------------------------"<<std::endl;
        trajectory = lane_change_trajectory(lane,state,sensor_fusion,car_info,
                                  previous_path_x,previous_path_y,end_path_s,end_path_d,
                                  map_waypoints_s,map_waypoints_x,
                                  map_waypoints_y,ref_val,
                                  intend_lane_info,final_lane_info);
    }else if(state == "PLCL" || state == "PLCR"){
        std::cout<<"------------------line 472: "<<state<<"------------------------"<<std::endl;
        trajectory = prep_lane_change_trajectory(lane,state,sensor_fusion,car_info,
                                  previous_path_x,previous_path_y,end_path_s,end_path_d,
                                  map_waypoints_s,map_waypoints_x,map_waypoints_y,ref_val,intend_lane_info,final_lane_info);
    }

    return trajectory;
}
float velocity_cost(float intend_lane_val,float final_lane_val){
    return 1.0-(intend_lane_val+final_lane_val)/(2*49.0);
}
float inlane_cost(double intend_lane,double final_lane){
    if(intend_lane >= 0 && intend_lane <= 2 && final_lane >=0 && final_lane <= 2)
        return 0;
    else
        return 1;
}
float calculate_cost(vector<float> intend_lane_info,vector<float> final_lane_info){
    float intend_lane_val = intend_lane_info[0];
    float intend_lane = intend_lane_info[1];
    float final_lane_val = final_lane_info[0];
    float final_lane = final_lane_info[1];
    float cost = 0.3*velocity_cost(intend_lane_val,final_lane_val)+
        0.7*inlane_cost(intend_lane,final_lane);
    return cost;
}
int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  int lane = 0;

  //double ref_val = 49;
  double ref_val = 0;
  string cur_state = "KL";

  h.onMessage([&cur_state,&ref_val,&lane,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

            vector<double> car_info{car_x,car_y,car_s,car_d,car_yaw,car_speed,end_path_s,end_path_d};
            lane = end_path_d/4;

            vector<double> next_x_vals;
            vector<double> next_y_vals;


            cout<<"last_state:"<<cur_state<<std::endl;
            vector<string> next_state = successor_states(cur_state,lane);

            double min_cost = -1;
            vector<point> best_trajectory;
            //string cur_state = "KL";
            double best_ref_val;
            for(auto state:next_state){
                double cur_ref_val = ref_val;
                vector<float> intend_lane_info;
                vector<float> final_lane_info;
                vector<point> trajectory = generate_trajectory(lane,state,
                        sensor_fusion,car_info,previous_path_x,
                        previous_path_y,end_path_s,end_path_d,
                        map_waypoints_s,map_waypoints_x,
                        map_waypoints_y,cur_ref_val,intend_lane_info,final_lane_info);
                float cost = calculate_cost(intend_lane_info,final_lane_info);

                std::cout<<"cur_state:"<<state<<std::endl;
                std::cout<<"ref_val"<<cur_ref_val<<std::endl;
                std::cout<<"intend_lane_val"<<intend_lane_info[0]<<std::endl;
                std::cout<<"intend_lane"<<intend_lane_info[1]<<std::endl;
                std::cout<<"final_lane_val"<<final_lane_info[0]<<std::endl;
                std::cout<<"final_lane"<<final_lane_info[1]<<std::endl;
                std::cout<<"cost:"<<cost<<std::endl;

                if(min_cost < 0 || cost < min_cost){
                    min_cost = cost;
                    best_trajectory = trajectory;
                    cur_state = state;
                    best_ref_val = cur_ref_val;
                }
            }
            ref_val = best_ref_val;
            cout<<"final_state:"<<cur_state<<std::endl;
            for(auto pose:best_trajectory){
                double x_point = pose.x;
                double y_point = pose.y;
                next_x_vals.push_back(x_point);
                next_y_vals.push_back(y_point);
            }


            //std::cout<<next_x_vals.size()<<std::endl;
            //for(int i=0;i<next_x_vals.size();i++){
                //std::cout<<next_x_vals[i]<<std::endl;
                //std::cout<<next_y_vals[i]<<std::endl;
            //}



            //Circle
            //
            //double pos_x;
            //double pos_y;
            //double angle;
            //int path_size = previous_path_x.size();
            //for(int i=0;i<path_size;i++){
                //next_x_vals.push_back(previous_path_x[i]);
                //next_y_vals.push_back(previous_path_y[i]);
            //}
            //if(path_size > 0){
                //pos_x = previous_path_x[path_size-1];
                //pos_y = previous_path_y[path_size-1];

                //double pos_x2 = previous_path_x[path_size-2];
                //double pos_y2 = previous_path_y[path_size-2];
                //angle = atan2(pos_y - pos_y2,pos_x - pos_x2);
            //}else if(path_size == 0){
                //pos_x = car_x;
                //pos_y = car_y;
                //angle = deg2rad(car_yaw);
            //}
            //double dist = 0.5;
            //for(int i=0;i<50-path_size;i++){
                //next_x_vals.push_back(pos_x + dist*cos(angle+(i+1)*pi()/100));
                //next_y_vals.push_back(pos_y + dist*sin(angle+(i+1)*pi()/100));

                //pos_x += dist*cos(angle+(i+1)*pi()/100);
                //pos_y += dist*sin(angle+(i+1)*pi()/100);
            //}

            //Follow the lane line
            //double dist = 0.5;
            //for(int i=0;i<50;i++){
                //double next_s = car_s + (i+1)*dist;
                //double next_d = 6;
                //vector<double> xy = getXY(next_s,next_d,map_waypoints_s,map_waypoints_x,map_waypoints_y);
                //next_x_vals.push_back(xy[0]);
                //next_y_vals.push_back(xy[1]);
                //std::cout<<xy[0]<<" "<<xy[1]<<std::endl;
            //}


          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	json msgJson;
          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
