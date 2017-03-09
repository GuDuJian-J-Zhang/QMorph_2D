//***********************************************************
// File:     USE_OpenGL.cpp
// Author:   Zhang Jun
// Date:     2016-05-12
// Email:    zhangjun_dg@mail.dlut.edu.cn
//***********************************************************

#include <Windows.h>
//#include <gl/glut.h>
#include <cmath>

#include "USE_OpenGL.h"

double Show_Mesh::_view_point_x = 0.0;
double Show_Mesh::_view_point_y = 0.0;

int Show_Mesh::_mouseX = 0;
int Show_Mesh::_mouseY = 0;
int Show_Mesh::_button = 0;
int Show_Mesh::_state = 0;
bool Show_Mesh::_mouseLeftDown = false;
bool Show_Mesh::_mouseRightDown = false;
int Show_Mesh::_currentH = 0;
int Show_Mesh::_currentW = 0;

using namespace std;

void Show_Mesh::ChangeSize(GLsizei w, GLsizei h)
{
	if (h == 0) h = 1;

	_currentW = w;
	_currentH = h;

	glViewport(50, -50, w, h);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	if (w <= h)
		glOrtho(-1.0f, 1.0f, -1.0f * h / w, 1.0f * h / w, 1.0f, -1.0f);
	else
		glOrtho(-1.0f * w / h, 1.0f * w / h, -1.0f, 1.0f, 1.0f, -1.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void Show_Mesh::SetupRC()
{
	glutReshapeWindow(1500,800);
	glutPositionWindow(400,200);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
}

void Show_Mesh::keyboard(unsigned char key ,int x,int y)
{
	if(key == 27)exit(0);
	else if(key == 'e')
		glutFullScreen();
	else if(key == 'd')
	{
		glutReshapeWindow(800,800);
		glutPositionWindow(400,200);
	}
}

void Show_Mesh::MouseFun(int button, int state, int x, int y)
{
	_mouseX = x;
	_mouseY = y;

	if (button == GLUT_LEFT_BUTTON)
	{
		if (state == GLUT_DOWN)
		{
			_mouseLeftDown = true;
		}			
		else
		{
			_mouseLeftDown = false;
		}
	}
	else if (button == GLUT_RIGHT_BUTTON)
	{
		if (state == GLUT_DOWN)
		{
			_mouseRightDown = true;
		}			
		else
		{
			_mouseRightDown = false;
		}
	}
}

void Show_Mesh::MouseMotion(int x, int y)
{
	if (_mouseLeftDown)//< ×ó¼ü
	{
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		if (_currentW <= _currentH)
			glOrtho(-1.0f, 1.0f, -1.0f * _currentH / _currentW, 1.0f * _currentH / _currentW, 1.0f, -1.0f);
		else
			glOrtho(-1.0f * _currentW / _currentH, 1.0f * _currentW / _currentH, -1.0f, 1.0f, 1.0f, -1.0f);

		glTranslatef((x - _mouseX)*0.005, (_mouseY - y)*0.005, 0);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glutPostRedisplay();
	}
	else if (_mouseRightDown) ///< ÓÒ¼ü
	{
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		if (_currentW <= _currentH)
			glOrtho(-1.0f, 1.0f, -1.0f * _currentH / _currentW, 1.0f * _currentH / _currentW, 1.0f, -1.0f);
		else
			glOrtho(-1.0f * _currentW / _currentH, 1.0f * _currentW / _currentH, -1.0f, 1.0f, 1.0f, -1.0f);
		glScalef(std::abs((_mouseX - x)), std::abs((_mouseX - x)), 1.0);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glutPostRedisplay();
	}
}

void Show_Mesh::show(void (*func)(void))
{
	glutInitDisplayMode(GLUT_SINGLE|GLUT_RGB);

	glutCreateWindow("Quad Mesh Generation (C) ZhangJun 2016");

	glutDisplayFunc(func);

	glutReshapeFunc(Show_Mesh::ChangeSize);

	glutKeyboardFunc(Show_Mesh::keyboard);

	glutMouseFunc(Show_Mesh::MouseFun);

	glutMotionFunc(Show_Mesh::MouseMotion);
	SetupRC();

	glutMainLoop();
}