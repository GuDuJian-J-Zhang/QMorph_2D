//**************************************************************************
// File:     foundation.h
// Author:   Zhang Jun
// Date:     2016-05-05
// Email:    zhangjun_dg@mail.dlut.edu.cn
//**************************************************************************

#ifndef _FOUNDATION_H
#define _FOUNDATION_H
#include "toolsdefine.h"
/******************************************************************************
*                                 Marker                                      *
******************************************************************************/

class Marker
{
public:
    Marker():_index(-1){}
public:
    void      setDead(bool dead=true) const
    {
        dead?_mark.set(0):_mark.reset(0);
    }
    bool      isDead()const
    {
        return _mark.test(0);
    }

	/** \brief 检查元素是否被做过某种标记
	*
	    \return true  元素已被标记
		\return false 元素未被标记
	*/
    bool      testMark(int pos) const 
    {
        return _mark.test(pos);
    }

	/** \brief 给元素做某种标记
	*
	    即将相应的标记位置1
	*/
    void      setMark(int pos) const
    {
        _mark.set(pos);
    }

	/** \brief 清除对元素的某种标记
	*
	    即将相应的标记位置0
	*/
    void      resetMark(int pos) const
    {
        _mark.reset(pos);
    }
    void      clearMark() const
    {
        _mark.reset();
    }
    void      setIndex(int idx)
    {
        _index = idx;
    }
    int       getIndex()const ///< 返回元素的索引号
    {
        return _index;
    }
private:
    mutable std::bitset<16>      _mark;  //元素标记 _maek有16个标记位（0~15），初始值均为0
    int                     _index; //元素编号
};

#endif
