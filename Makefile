raycast:
	g++ -g -Wall -Wno-deprecated ray.cpp glm.cpp udray.cpp -o RayCast -lglut -lGLU -lGL

test:
	g++ test.cpp ray-test.cpp glm.cpp udray.cpp -o Test -lglut -lGLU -lGL
