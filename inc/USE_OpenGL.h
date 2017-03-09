//***********************************************************
// File:     USE_OpenGL.h
// Author:   Zhang Jun
// Date:     2016-05-12
// Email:    zhangjun_dg@mail.dlut.edu.cn
//***********************************************************
#include <glew.h>
#include <freeglut.h>

class Show_Mesh
{
public:
	Show_Mesh() {}
	~Show_Mesh() {}
	void show(void (*func)(void));
	void set_view_point(double x, double y) {_view_point_x = x; _view_point_y = y;}

private:
	static void ChangeSize(GLsizei w, GLsizei h);
	static void SetupRC(void);
	static void keyboard(unsigned char key ,int x,int y);
	static void MouseFun(int button, int state, int x, int y);
	static void MouseMotion(int x, int y);
	//void SpecialKeys(int key, int x, int y);

private:
	static double _view_point_x;
	static double _view_point_y;
	static int    _mouseX;
	static int    _mouseY;
	static int    _button;
	static int    _state;
	static bool   _mouseLeftDown;
	static bool   _mouseRightDown;
	static int    _currentH;
	static int    _currentW;

};