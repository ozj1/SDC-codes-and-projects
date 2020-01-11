#include <uWS/uWS.h>
#include <iostream>
#include "json.hpp"
#include "PID.h"
#include <math.h>
#include <deque>

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() {
  return M_PI;
}
double deg2rad(double x) {
  return x * pi() / 180;
}
double rad2deg(double x) {
  return x * 180 / pi();
}

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
std::string hasData(std::string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_last_of("]");
  if (found_null != std::string::npos) {
    return "";
  } else if (b1 != std::string::npos && b2 != std::string::npos) {
    return s.substr(b1, b2 - b1 + 1);
  }
  return "";
}

int main() {
  uWS::Hub h;

  // steering pid
  PID pid;
  // TODO: Initialize the pid variable.
  double Kp = 0.1;  // proportional coefficient
  double Ki = 0.003;  // integral coefficient
  double Kd = 3.0;  // differential coefficient
  pid.Init(Kp, Ki, Kd);

  std::deque<double> steering_history;

  // throttle pid
  PID pidt;
  pidt.Init(2, 0.001, 4);

  h.onMessage(
      [&pid,&pidt,&steering_history](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length, uWS::OpCode opCode) {
        // "42" at the start of the message means there's a websocket message event.
        // The 4 signifies a websocket message
        // The 2 signifies a websocket event
        if (length && length > 2 && data[0] == '4' && data[1] == '2')
        {
          auto s = hasData(std::string(data));
          if (s != "") {
            auto j = json::parse(s);
            std::string event = j[0].get<std::string>();
            if (event == "telemetry") {
              // j[1] is the data JSON object
              double cte = std::stod(j[1]["cte"].get<std::string>());
              double speed = std::stod(j[1]["speed"].get<std::string>());
              double angle = std::stod(j[1]["steering_angle"].get<std::string>());
              double steer_value;
              /*
               * TODO: Calculate steering value here, remember the steering value is
               * [-1, 1].
               * NOTE: Feel free to play around with the throttle and speed. Maybe use
               * another PID controller to control the speed!
               */


              // steering history
              steering_history.push_front(angle);

              // only keep so many
              if (steering_history.size() > 20)
                steering_history.pop_back();

              double avg_angle = std::accumulate(steering_history.begin(), steering_history.end(), 0.0) / steering_history.size();
              double abs_avg_angle = fabs(avg_angle);

              // setup desired speed
              double desired_speed= 50;

              // get some speed up first
              if (speed > 15) {
                if (abs_avg_angle < 1)
                  desired_speed = 100;
                else if (abs_avg_angle < 1.5)
                  desired_speed = 90;
                else if (abs_avg_angle < 3)
                  desired_speed = 65;
                else if (abs_avg_angle < 4)
                  desired_speed = 55;
                else if (abs_avg_angle < 5)
                  desired_speed = 45;
                else if (abs_avg_angle < 6)
                  desired_speed = 35;
                else if (abs_avg_angle < 8)
                  desired_speed = 25;
                else if (abs_avg_angle > 10)
                  desired_speed = 15; // slow down
              }

              // adjust desired speed also for CTE
              if (speed > 25) {
                double acte = fabs(cte);
                if (acte < .1)
                  desired_speed += 20;
                else if (acte < .2  )
                  desired_speed += 10;
                else if (acte < .3)
                  desired_speed -= 5;
                else if (acte < .5 )
                  desired_speed -= 15;
                else if (acte < .7)
                  desired_speed -= 25;
                else if (acte > .8)
                  desired_speed = -10; // break
                else
                  desired_speed = 5;// slow down
              }
              pidt.UpdateError(speed-desired_speed);
              double throttle = pidt.TotalError();


              // update steering
//              if (speed > 55) {
//                pid.Kp = 0.1;
//                pid.Ki = 0.004;
//                pid.Kd = 6;
//              } else {
//                pid.Kp = 0.1;
//                pid.Ki = 0.004;
//                pid.Kd = 3;
//              }
              pid.UpdateError(cte);
              steer_value = pid.TotalError();
              // DEBUG
              std::cout << "CTE: " << cte << " Steering Value: " << steer_value << " Angle: " << angle << " Avg Angle: " << avg_angle
              << " Speed: " << speed << std::endl;

              json msgJson;
              msgJson["steering_angle"] = steer_value;
              //msgJson["throttle"] = 0.3;
              msgJson["throttle"] = throttle;

              auto msg = "42[\"steer\"," + msgJson.dump() + "]";
              std::cout << msg << std::endl;
              ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
            }
          } else {
            // Manual driving
            std::string msg = "42[\"manual\",{}]";
            ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          }
        }
      });

  // We don't need this since we're not using HTTP but if it's removed the program
  // doesn't compile :-(
  h.onHttpRequest(
      [](uWS::HttpResponse *res, uWS::HttpRequest req, char *data, size_t, size_t) {
        const std::string s = "<h1>Hello world!</h1>";
        if (req.getUrl().valueLength == 1)
        {
          res->end(s.data(), s.length());
        }
        else
        {
          // i guess this should be done more gracefully?
          res->end(nullptr, 0);
        }
      });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection(
      [&h](uWS::WebSocket<uWS::SERVER> ws, int code, char *message, size_t length) {
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
