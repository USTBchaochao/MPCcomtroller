import MPCcontroller as MPC
import numpy
####模拟环境y=2直线
N = 100
Route = []
Ts = 0.05
for i in range(N):
    Route.append([i*Ts , 2 , 0])
# print(numpy.shape(Route)[0])

##模拟车辆
State = numpy.array(([0],[1],[2]))#初始状态【x，y，航向角】
Speed = numpy.array(([3],[4]))#初始速度【速度，前轮转角】
L = 1
##初始化
Np = 4#预测时域
Nc = 3#控制时域
ROW = numpy.array(([1]))
dumin = -10##速度增量限制
dumax = 10

dsmin = -0.5#转向增量限制
dsmax = 0.5
dmin = 0##松弛因子限制
dmax = 100

dspeed , dsteer = MPC.MPCController(State , Speed , Ts , Np , Nc , dumax , dumin , dsmax , dsmin , dmax , dmin , ROW , L , Route)
speed += dspeed
steer +=dsteer
print(commandspeed , commandsteer)
# print(result['x'][0])
