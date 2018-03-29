import os 

def readToMas(path):
	mas = []
	tmp = os.listdir(os.getcwd()+ path)
	for fileName in tmp:
		file = open(os.getcwd()+ path + "/" +fileName, "r")
		string = file.read().split()
		mas.append((float(string[2]), float(string[11]), float(string[15]), float(string[19])))
		file.close()
	mas = sorted(mas, key=lambda x: x[0])
	return mas;

def calc(mas):
	time = []
	accel = []
	eff = []
	for i in range(len(mas)):
		time.append(mas[i][1] + mas[i][2] + mas[i][3])
	t = mas[0][1] + mas[0][2] + mas[0][3]
	for i in range(1, len(mas)):
		s = t / (float(mas[i][1]) + float(mas[i][2]) + float(mas[i][3]))
		accel.append(s)
		eff.append(s / float(mas[i][0]))
	return (time, accel, eff)

def tableTime(xAxis, mas):
	for i in range(len(mas)):
		print(str(int(mas[i][0])) + " & " + str(mas[i][1]) + " & " + str(mas[i][2]) + " & " + str(mas[i][3]) + " \\\\")

def tableAccel(xAxis, accel):
	for i in range(len(accel)):
		print(str(xAxis[i]) + " & " + str(accel[i]) + " \\\\")

def tableEff(xAxis, eff):
	for i in range(len(eff)):
		print(str(xAxis[i]) + " & " + str(eff[i]) + " \\\\")

def graphTime(mas1024, mas2048, mas3072):
	for i in range(0,len(mas1024)):
		print(str(i+1) + " " + str(mas1024[i]) + " " + str(mas2048[i]) + " " + str(mas3072[i]))

mas1024 = readToMas("/x1024")
mas2048 = readToMas("/x2048")
mas3072 = readToMas("/x3072")

(time1024, accel1024, eff1024) = calc(mas1024)
(time2048, accel2048, eff2048) = calc(mas2048)
(time3072, accel3072, eff3072) = calc(mas3072)

xAxis = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]

# tableEff(xAxis, eff3072)

# tableEff(xAxis, eff1024)
graphTime(eff1024, eff2048, eff3072)
