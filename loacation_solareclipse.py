#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 20:54:21 2018

@author: wpf
"""
from novas import compat as novas
from novas.compat  import eph_manager
from novas  import constants
import math
#import datetime
#时间初始值设置
year=2013
month=11
day=3
hour=10
#读取计算机时间
'''
nowtime=datetime.datetime.now().strftime('%Y-%m-%d %H%M%S')
year=int(nowtime[0:4])
month=int(nowtime[5:7])
day=int(nowtime[8:10])-1
hour=int(nowtime[11:13])+int(nowtime[13:15])/60+int(nowtime[15:17])/3600-8
'''
#儒略日转成北京时间
def jdutc2bt(jd_utc):
    t=novas.cal_date(jd_utc+8/24 )
    a = int(t[3])
    b = int((t[3] - a) * 60)
    c = (t[3] - a - b / 60.0)*3600.0
    return (t,a,b,c) 
#地理位置
    #佘山
latitude=31.084451#纬度
longitude=121.165869#经度
height=12.2#海拔
    #老家
'''latitude=34.173863#纬度
longitude=115.927278#经度
height=42#海拔'''

tempetature=25.0#温度
pressure=1010.0#气压
#距离以AU为单位
#re = (constants.ERAD + 65000.0 )/ constants.AU#地球半径，考虑地球大气层的影响
re = (constants.ERAD) / constants.AU#地球半径，不考虑地球大气层的影响
rs = 696000000.0 / constants.AU#太阳半径
rm = 1738000.0 / constants.AU#月亮半径
ta = constants.AU / constants.C#从太阳到地球的光行时
#判断地球上某一地点是否是白天
def day_time(ETT,E_O4):
    vpd=0
    ETTlength=0
    E_O4length=0
    for i in range(3):
        vpd+=ETT[i]*E_O4[i]
        ETTlength+=ETT[i]*ETT[i]
        E_O4length+=E_O4[i]*E_O4[i]
    ETTlength=math.sqrt(ETTlength)
    E_O4length=math.sqrt(E_O4length)
    theta=math.acos(vpd/(ETTlength*E_O4length))
    if theta>(math.pi/2-math.asin(rs-re)):
        return True
    else:
        return False
dtheta=0.0
leap_second=37#闰秒
tt_tai=32.184#tt和国际原子时tai之差
#EOP参数
ut1_utc=0.06723
x_pole = -0.002
y_pole = +0.529
#tt-ut1
delta_t=tt_tai+leap_second-ut1_utc
#构造不同的儒略日时间变量
jd_utc=novas.julian_date(year,month,day,hour)#utc
jd_ut1=jd_utc+ut1_utc/86400.0#ut1
jd_tt=jd_utc+(leap_second+tt_tai)/86400.0#tt
jd=(jd_tt,0.0)
jd0=(jd_tt-ta/86400,0.0)
#打开de历表
jd_s,jd_e,num=eph_manager.ephem_open()
 #太阳和地球构造
sun=novas.make_object(0,10,'sun',None)
moon=novas.make_object(0,11,'moon',None)
earth=novas.make_object(0,3,'earth',None)
 #位置构造
location1=novas.make_on_surface(latitude,longitude,height,25.0,1013)
location=novas.make_observer_on_surface(latitude,longitude,height,25.0,1013)
#矢量和夹角初始化
O1_S1=[0.0,0.0,0.0]#月球半影锥点到太阳质心矢量坐标
O2_S1=[0.0,0.0,0.0]#月球全影锥点到太阳质心矢量坐标
O1_E=[0.0,0.0,0.0]#月球半影锥点到地心距离矢量坐标
O1_T=[0.0,0.0,0.0]#月球半影锥点到地球某一点距离矢量坐标
O2_E=[0.0,0.0,0.0]#月球全影锥点到地心距离矢量坐标
O2_T=[0.0,0.0,0.0]#月球全影锥点到地球某一点距离矢量坐标
O1_M=[0.0,0.0,0.0]#月球半影锥点到月心距离矢量坐标
O2_M=[0.0,0.0,0.0]#月球全影锥点到月心距离矢量坐标
E_O4=[0.0,0.0,0.0]
#经纬度转为ITRS下坐标
E_T=[re*math.cos(latitude*constants.DEG2RAD)*math.cos(longitude*constants.DEG2RAD),\
	 re*math.cos(latitude*constants.DEG2RAD)*math.sin(longitude*constants.DEG2RAD),\
	 re*math.sin(latitude*constants.DEG2RAD)\
     ]
#ITRS->GCRS
E_TT=novas.ter2cel(jd_ut1, 0.0, delta_t, x_pole, y_pole,E_T, 1, 0,0)#
thetaE_O1=0.0#地球半径在半影锥点O1角距
thetaM_O1=0.0#月球半径在半影锥点O1角距
thetaE_O2=0.0#地球半径在全影锥点O2角距
thetaM_O2=0.0#月球半径在全影锥点O2角距
thetaEM_O1=10.0#地心和月心在月球半影锥点为原点的夹角
thetaTM_O1=10.0#地球某一点和月心在月球半影锥点为原点的夹角
thetaEM_O2=1.0#地心和月球在月心全影锥点为原点的夹角
thetaTM_O2=1.0#地球某一点和月球在月心全影锥点为原点的夹角
O2_Elength=1#月球全影锥点可能在地球内部，需要判断，先赋初值足够大
O2E_re=1#月球全影锥点到地心距离和地球半径之差
#theta=constants.TWOPI
flag1=0#日偏食开始结束标志，用来结束循环
flag2=0#日全食发生标记
flag3=0#日环食标记
flag4=0#食甚标记
flag_a=0#日环食时月球全影锥点进出地球标记
#某一点
flag11=0#日偏食开始结束标志，用来结束循环
flag22=0#日全食发生标记
flag33=0#日环食标记
flag44=0#食甚标记
while True:
    #partial
    dthetaEM_O1=thetaEM_O1-thetaE_O1-thetaM_O1
    dthetaTM_O1=thetaTM_O1-thetaM_O1
    #total
    dthetaEM_O2=thetaEM_O2-thetaE_O2-thetaM_O2
    dthetaTM_O2=thetaTM_O2-thetaM_O2
    #annular
    dthetaEM_O2A=thetaEM_O2+thetaE_O2+thetaM_O2-math.pi
    dthetaTM_O2A=thetaTM_O2+thetaM_O2-math.pi
    O2E_re=math.fabs(O2_Elength-re)#在地球内部设其值为1
    thetaEM=thetaEM_O1#上一步的地月夹角，食甚时夹角最小
    thetaTM=thetaTM_O1#上一步的地月夹角，食甚时夹角最小
    jd_utc+=0.1/86400#每次增加0.1秒
    jd_ut1=jd_utc+ut1_utc/86400.0#ut1
    jd_tt=jd_utc+(leap_second+tt_tai)/86400.0#tt
    jd=(jd_tt,0.0)#tt代替tdb，差别不大
    jd0=(jd_tt-ta/86400,0.0)#太阳发出光时的时间
    pos_earth0=novas.ephemeris(jd0,earth)#太阳发出光时的地球icrs坐标
    pos_earth=novas.ephemeris(jd,earth)#太阳到达地球时的地球icrs坐标
    pos_moon0=novas.ephemeris(jd0,moon)#太阳发出光时的月球icrs坐标
    pos_moon=novas.ephemeris(jd,moon)#太阳到达地球时的月球icrs坐标
    vpd_O1EM=0.0#O1E和O1M矢量积
    vpd_O1TM=0.0#O1T和O1M矢量积
    vpd_O2EM=0.0#O2E和O2M矢量积
    vpd_O2TM=0.0#O2T和O2M矢量积
    O1_Elength=0.0#月球半影锥点到地心距离
    O1_Tlength=0.0#月球半影锥点到地球某一点距离
    O1_Mlength=0.0#月球半影锥点到月心距离
    O2_Elength=0.0#月球全影锥点到地心距离
    O2_Tlength=0.0#月球全影锥点到地球某一点距离
    O2_Mlength=0.0#月球全影锥点到月心距离
    #ITRS->GCRS
    E_TT=novas.ter2cel(jd_ut1, 0.0, delta_t, x_pole, y_pole,E_T, 1, 0,0)#
    for i in range(3):
        O1_S1[i] = rm / (rm+rs) * pos_moon0[0][i]-pos_moon[0][i]
        O2_S1[i] = rm/  (rm-rs) * pos_moon0[0][i]-pos_moon[0][i] 
        O1_M[i] = O1_S1[i] + pos_moon[0][i]
        O2_M[i] = O2_S1[i] + pos_moon[0][i]
        O1_E[i] = O1_S1[i] + pos_earth[0][i]
        O1_T[i] = O1_E[i]+E_TT[i]#月球半影锥点到地球上某一点的矢量
        O2_E[i] = O2_S1[i] + pos_earth[0][i]
        O2_T[i] = O2_E[i]+E_TT[i]#月球全影锥点到地球上某一点的矢量
        E_O4[i] = re/(rs-re)*pos_earth0[0][i]#用于判断某一点是否处于白天
        O1_Elength+=O1_E[i]*O1_E[i]
        O1_Tlength+=O1_T[i]*O1_T[i]#月球半影锥点到地球上某一点的矢量长度
        O1_Mlength+=O1_M[i]*O1_M[i]
        O2_Elength+=O2_E[i]*O2_E[i]
        O2_Tlength+=O2_T[i]*O2_T[i]#月球全影锥点到地球上某一点的矢量长度
        O2_Mlength+=O2_M[i]*O2_M[i]
        vpd_O1EM+=O1_E[i]*O1_M[i]
        vpd_O1TM+=O1_T[i]*O1_M[i]
        vpd_O2EM+=O2_E[i]*O2_M[i]
        vpd_O2TM+=O2_T[i]*O2_M[i]
    #矢量的长度
    O1_Elength=math.sqrt(O1_Elength)#O1_E长度
    O1_Tlength=math.sqrt(O1_Tlength)#O1_T长度
    O1_Mlength=math.sqrt(O1_Mlength)#O1_M长度
    O2_Elength=math.sqrt(O2_Elength)#O2_E长度
    O2_Tlength=math.sqrt(O2_Tlength)#O2_T长度
    O2_Mlength=math.sqrt(O2_Mlength)#O2_M长度
    #地球和月球在O1,O2点为原点的半径的角距
    thetaE_O1=math.asin(re/O1_Elength)
    thetaM_O1=math.asin(rm/O1_Mlength)
    if O2_Elength<re:#全影锥点在地球内部
        thetaE_O2=math.pi/2#在地球内部设其值为1
    else:
        thetaE_O2=math.asin(re/O2_Elength)
    thetaM_O2=math.asin(rm/O2_Mlength)
    #矢量O1E,O1M的夹角和矢量O2E,O2M的夹角
    thetaEM_O1=math.acos(vpd_O1EM/(O1_Elength*O1_Mlength))
    thetaTM_O1=math.acos(vpd_O1TM/(O1_Tlength*O1_Mlength))
    thetaEM_O2=math.acos(vpd_O2EM/(O2_Elength*O2_Mlength))
    thetaTM_O2=math.acos(vpd_O2TM/(O2_Tlength*O2_Mlength))
    #日偏食
    if dthetaEM_O1*(thetaEM_O1-thetaE_O1-thetaM_O1)<0: #（1）flag1==flag4==0：日偏食开始之前
        if flag1==0:                                                             #（2）flag1==flag4==1：食甚之后，日偏食结束之前
            t1=jdutc2bt(jd_utc)#开始
            flag1+=1
        elif flag1==1:
            t2=jdutc2bt(jd_utc)#结束
            break#日偏食结束，日食结束
    if day_time(E_TT,E_O4):#先判断是否处于白天
        if dthetaTM_O1*(thetaTM_O1-thetaM_O1)<0: #（1）flag1==flag4==0：日偏食开始之前
            if flag11==0:                                                             #（2）flag1==flag4==1：食甚之后，日偏食结束之前
                t11=jdutc2bt(jd_utc)#开始
                flag11+=1
            elif flag11==1:
                t22=jdutc2bt(jd_utc)#结束
                flag11+=1
     #日环食
    if dthetaEM_O2A*(thetaEM_O2+thetaE_O2+thetaM_O2-math.pi)*O2E_re*(O2_Elength-re)<0 and flag_a==0:#两个判断条件不能同时满足 ，防止后者引发前者的发生
        if flag3==0:
            t5=jdutc2bt(jd_utc)#开始
            flag3+=1
        elif flag3>0:
            t6=jdutc2bt(jd_utc)#结束
            flag3+=1
        if O2E_re*(O2_Elength-re)<0:
            flag_a+=1
    if day_time(E_TT,E_O4):
        if dthetaTM_O2A*(thetaTM_O2+thetaM_O2-math.pi)<0:#和日全食类似
            if flag33==0:
                t55=jdutc2bt(jd_utc)#开始
                flag33+=1
            elif flag33>0:
                t66=jdutc2bt(jd_utc)#结束
                flag33+=1
    #日全食
    if dthetaEM_O2*(thetaEM_O2-thetaE_O2-thetaM_O2)<0 or (O2E_re*(O2_Elength-re)<0):
        if flag2==0:
            t3=jdutc2bt(jd_utc)#开始
            flag2+=1
        elif flag2>0:
            t4=jdutc2bt(jd_utc)#结束
            flag2+=1
    if day_time(E_TT,E_O4):
        if dthetaTM_O2*(thetaTM_O2-thetaM_O2)<0:#和日偏食类似
            if flag22==0:
                print("fuck")
                t33=jdutc2bt(jd_utc)#开始
                flag22+=1
            elif flag22>0:
                t44=jdutc2bt(jd_utc)#结束
    #食甚
    if thetaEM<thetaEM_O1 and flag4==0:#食甚按照半影影锥锥点来算，因为全影影锥和半影影锥的大小变化、同时达到最值的时间大小不一定相同
        t7=jdutc2bt(jd_utc)#食甚
        flag4+=1
    if day_time(E_TT,E_O4):
        if thetaTM<thetaTM_O1 and flag44==0:#食甚按照半影影锥锥点来算，因为全影影锥和半影影锥的大小变化、同时达到最值的时间大小不一定相同
            t77=jdutc2bt(jd_utc)#食甚
            flag44+=1
#换算成北京时间输出
#日食总体情况
print("{0}-{1}-{2}号的日食：".format(t1[0][0], t1[0][1], t1[0][2]))
if flag2==0 and flag3==0:
    print("日食类型：日偏食")
if flag2!=0 and flag3==0:
    print("日食类型：日全食")
if flag2==0 and flag3!=0:
    print("日食类型：日环食")
if flag2>1 and flag3>1:#日全食和日环食都有有始有终，若是只有开始，则可能是情况4
    print("日食类型：全环食")
print("日偏食开始：{0:02d}:{1:02d}:{2:4.2f}".format(t1[1], t1[2], t1[3]))
if flag3>1:
    print("日环食开始：{0:02d}:{1:02d}:{2:4.2f}".format(t5[1], t5[2], t5[3]))
if flag2!=0:
    print("日全食开始：{0:02d}:{1:02d}:{2:4.2f}".format(t3[1], t3[2], t3[3]))
print("食     甚：{0:02d}:{1:02d}:{2:4.2f}".format(t7[1], t7[2], t7[3]))
if flag2!=0:
    print("日全食结束：{0:02d}:{1:02d}:{2:4.2f}".format(t4[1], t4[2], t4[3]))
if flag3>1:
    print("日环食结束：{0:02d}:{1:02d}:{2:4.2f}".format(t6[1], t6[2], t6[3]))
print("日偏食结束：{0:02d}:{1:02d}:{2:4.2f}".format(t2[1], t2[2], t2[3]))
#某一地点的日食情况
print("经纬度为（{0},{1}）的地方的日食：".format(longitude,latitude))
 #日食类型判断
if flag11==0 and flag22==0 and flag33==0:
   print("没有日食发生")
if flag11!=0 and flag22==0 and flag33==0:
    print("日食类型：日偏食")
if flag22!=0 and flag33==0:
    print("日食类型：日全食")#
if flag22==0 and flag33!=0:
    print("日食类型：日环食")
if flag22!=0 and flag33!=0:
    print("日食类型：全环食")
 #具体时间输出
if flag11==1:#满足日偏食开始或者结束的条件只发生了一次，意味着日食还未结束就已经日落后或者日出时已经不是日偏食了
    if t11[1]>12:#日偏食发生时间在下午，意味着带食日落
        print("带食日落")#具体带什么食日落懒得判断了，日食过程如下
        print("日偏食开始：{0:02d}:{1:02d}:{2:4.2f}".format(t11[1], t11[2], t11[3]))
        if flag22>0:#意味着有日全食开始
            print("日全食开始：{0:02d}:{1:02d}:{2:4.2f}".format(t33[1], t33[2], t33[3]))
        if flag33>0:#意味着有日环食开始
            print("日环食开始：{0:02d}:{1:02d}:{2:4.2f}".format(t55[1], t55[2], t55[3])) 
        if flag44!=0:#意味着有食甚发生
            print("食     甚：{0:02d}:{1:02d}:{2:4.2f}".format(t77[1], t77[2], t77[3]))
        if flag22>1:#意味着有日全食结束
            print("日全食结束：{0:02d}:{1:02d}:{2:4.2f}".format(t44[1], t44[2], t44[3]))
        if flag33>1:#意味着有日环食结束
            print("日环食结束：{0:02d}:{1:02d}:{2:4.2f}".format(t66[1], t66[2], t66[3]))
    else:#日偏食发生时间在上午，意味着带食带食日出
        print("带食日出")#具体带什么食日出懒得判断了，日食过程如下
        if flag22>1:#意味着有日全食开始
            print("日全食开始：{0:02d}:{1:02d}:{2:4.2f}".format(t33[1], t33[2], t33[3]))
        if flag22>0:#意味着有日全食开始
            if flag22>1:#意味着有日全食结束，也意味着有食甚
                print("日全食开始：{0:02d}:{1:02d}:{2:4.2f}".format(t33[1], t33[2], t33[3]))
                print("食     甚：{0:02d}:{1:02d}:{2:4.2f}".format(t77[1], t77[2], t77[3]))
                print("日全食结束：{0:02d}:{1:02d}:{2:4.2f}".format(t44[1], t44[2], t44[3]))
            else:
                print("日全食结束：{0:02d}:{1:02d}:{2:4.2f}".format(t33[1], t33[2], t33[3]))
        if flag33>0:#意味着有日环食开始
            if flag33>1:#意味着有日环食结束，也意味着有食甚
                print("日环食开始：{0:02d}:{1:02d}:{2:4.2f}".format(t55[1], t55[2], t55[3])) 
                print("食     甚：{0:02d}:{1:02d}:{2:4.2f}".format(t77[1], t77[2], t77[3]))
                print("日环食结束：{0:02d}:{1:02d}:{2:4.2f}".format(t66[1], t66[2], t66[3]))
            else:
                print("日环食结束：{0:02d}:{1:02d}:{2:4.2f}".format(t55[1], t55[2], t55[3])) 
        print("日偏食结束：{0:02d}:{1:02d}:{2:4.2f}".format(t11[1], t11[2], t11[3]))
if flag11>1:#满足日偏食开始或者结束的条件只发生了两次，日食过程中没有日出日落
    print("日偏食开始：{0:02d}:{1:02d}:{2:4.2f}".format(t11[1], t11[2], t11[3]))
    if flag22==0 and flag33==0:
        print("食     甚：{0:02d}:{1:02d}:{2:4.2f}".format(t77[1], t77[2], t77[3]))
    if flag22!=0:#因为日食过程中没有日出日落，所以有始有终
        print("日全食开始：{0:02d}:{1:02d}:{2:4.2f}".format(t33[1], t33[2], t33[3]))
        print("食     甚：{0:02d}:{1:02d}:{2:4.2f}".format(t77[1], t77[2], t77[3]))
        print("日全食结束：{0:02d}:{1:02d}:{2:4.2f}".format(t44[1], t44[2], t44[3]))
    if flag33!=0:#因为日食过程中没有日出日落，所以有始有终
        print("日环食开始：{0:02d}:{1:02d}:{2:4.2f}".format(t55[1], t55[2], t55[3]))
        print("食     甚：{0:02d}:{1:02d}:{2:4.2f}".format(t77[1], t77[2], t77[3]))
        print("日环食结束：{0:02d}:{1:02d}:{2:4.2f}".format(t66[1], t66[2], t66[3]))
    print("日偏食结束：{0:02d}:{1:02d}:{2:4.2f}".format(t22[1], t22[2], t22[3]))
