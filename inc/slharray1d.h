//slharray1d.h
#ifndef __SLHARRAY1D_H_
#define __SLHARRAY1D_H_
#include <sldelaunay2d_std.h>
SL_MSH_NS_BEGIN
template <class T> 
class ObjectCopy
{
public:
    static void Memory_copy(T* pCopy, const T * pSource, int nCount);
    static void Memory_init(T * mem, int nCount,const T& val=T());
};
template <class T,int IS=1024,class R = ObjectCopy<T> >
class HArray1D
{
public:
    HArray1D(UInt sz=0);
    HArray1D(const T* mem,UInt sz);
    HArray1D(UInt sz,const T& val);
    virtual ~HArray1D(){clear();}

    HArray1D<T,IS,R>& operator=(const T& sc);

    UInt  size()const{return _size;}
    UInt  capacity()const{return _physize;}
    void  clear();
    bool  empty()const{return _size==0;}
    void  swap(HArray1D<T,IS,R>& sc);

    void            reset(const T& val);
    T&              operator[](UInt index){return at(index);}
    const T&        operator[](UInt index)const{return at(index);}
    T&              at(UInt index);
    const T&        at(UInt index)const;

    void            reserve(UInt sz,const T& val=T());
    void            resize(UInt sz,const T& val=T());
    void            push_back(const T&);
    void            pop_back();

    typedef Tuple<T,IS>      ATuple;
    typedef Array1D<ATuple*> ATupleHeads;
    typedef T*       pointer;
    typedef const T* const_pointer;
    typedef const T& const_reference;
    typedef T&       reference;
    struct iterator 
    {
        UInt  _whi,_offset;
        ATupleHeads* _data;
        iterator(UInt whi,UInt offset,ATupleHeads* data) 
            : _whi(whi), _offset(offset),_data(data) {}

        iterator() {}
        reference operator*() const 
        {
            return (*((*_data)[_whi]))[_offset]; 
        }
        pointer operator->() const { return &(operator*()); }
        iterator& operator++()
        {
            if(_offset>=IS)
            {
                _offset=0;
                _whi+=1;
            }else
            {
                _offset+=1;
            }
            return *this;
        }

        iterator operator++(int)
        {
            iterator tmp = *this;
            ++*this;
            return tmp;
        }

        bool operator==(const iterator& it) const
        { 
            return (_offset==it._offset &&  _whi==it._whi);
        }
        bool operator!=(const iterator& it) const
        {
            return (_offset!=it._offset ||  _whi!=it._whi);
        }
    };

    struct const_iterator 
    {
        UInt  _whi,_offset;
        ATupleHeads* _data;
        const_iterator(UInt whi,UInt offset,ATupleHeads* data) 
            : _whi(whi), _offset(offset),_data(data) {}

        const_iterator() {}
        const_reference operator*() const { return (*((*_data)[_whi]))[_offset];}
        const_pointer operator->() const { return &(operator*()); }
        const_iterator& operator++()
        {
            if(_offset>=IS)
            {
                _offset=0;
                _whi+=1;
            }else
            {
                _offset+=1;
            }
            return *this;
        }

        const_iterator operator++(int)
        {
            iterator tmp = *this;
            ++*this;
            return tmp;
        }

        bool operator==(const const_iterator& it) const
        { 
            return (_offset==it._offset &&  _whi==it._whi);
        }

        bool operator!=(const const_iterator& it) const
        {
            return (_offset!=it._offset ||  _whi!=it._whi);
        }
    };

    iterator        begin(){return iterator(0,0,&_data);}
    const_iterator  begin()const{return const_iterator(0,0,&_data);}
    iterator        end()
    {
        UInt pos=_size/IS;
        UInt off=_size%IS;
        return iterator(pos,off,&_data);
    }
    const_iterator  end()const
    {
        UInt pos=_size/IS;
        UInt off=_size%IS;
        return const_iterator(pos,off,&_data);
    }
    void  erase(iterator a,iterator b);
    void  erase(const_iterator a,const_iterator b);
private:
    HArray1D(const HArray1D<T,IS,R>& sc){}
    HArray1D<T,IS,R>& operator=(const HArray1D<T,IS,R>& sc){return *this;}
protected:
    ATupleHeads   _data;
    UInt          _size,_physize;
};

template <class T> 
void ObjectCopy<T>:: Memory_copy(T* pCopy, const T * pSource, int nCount)
{
    for(int i=0;i<nCount;i++)
    {
        pCopy[i]=pSource[i];
    }
}

template <class T> 
void ObjectCopy<T>::Memory_init(T * mem, int nCount,const T& val)
{
    for (int i=0;i<nCount;i++) 
    {
        mem[i]=val;
    }
}

template <class T,int IS,class R>
HArray1D<T,IS,R>::HArray1D(UInt sz)
:_size(0),_physize(0)
{
    resize(sz);
}



template <class T,int IS,class R>
HArray1D<T,IS,R>& HArray1D<T,IS,R>::operator=(const T& sc)
{
    for(UInt i=0;i<_size;i++)
    {
        this->at(i)=ac;
    }
}


template <class T,int IS,class R>
HArray1D<T,IS,R>::HArray1D(const T* mem,UInt sz)
:_size(0),_physize(0)
{
    resize(sz);
    for(UInt i=0;i<sz;i++)
    {
        this->at(i)=mem[i];
    }
}

template <class T,int IS,class R>
HArray1D<T,IS,R>::HArray1D(UInt sz,const T& val)
:_size(0),_physize(0)
{
    resize(sz,val);
}

template <class T,int IS,class R>
void HArray1D<T,IS,R>::clear()
{
    int i,num=_physize/IS;
    for(i=0;i<num;i++)
    {
        delete _data[i];
    }
    _physize=_size=0;
    _data.clear();
}

template <class T,int IS,class R>
void  HArray1D<T,IS,R>::swap(HArray1D<T,IS,R>& sc)
{
    Swap(this->_data,sc._data);
    Swap(_size,sc._physize);
    Swap(_physize,sc._physize);
}

template <class T,int IS,class R>
T&  HArray1D<T,IS,R>::at(UInt index)
{
    _MTAssert(index<_size);
    UInt pos=index/IS;
    UInt off=index%IS;
    return (*(_data[pos]))[off];
}

template <class T,int IS,class R>
const T& HArray1D<T,IS,R>::at(UInt index)const
{
    _MTAssert(index<_size);
    UInt pos=index/IS;
    UInt off=index%IS;
    return (*(_data[pos]))[off];
}

template <class T,int IS,class R>
void HArray1D<T,IS,R>::push_back(const T& val)
{
    if(_size>=_physize)
    {
        ATuple* anew=SL_NEW ATuple();
        R::Memory_init(anew->data,IS,val);
        _data.push_back(anew);
        _physize+=IS;
    }
    UInt pos=_size/IS;
    UInt off=_size%IS;
    (*(_data[pos]))[off]=val;
    _size++;
}

template <class T,int IS,class R>
void HArray1D<T,IS,R>::pop_back()
{
    _size-=1;
    UInt diff=_physize-_size;
    if(diff>IS)
    {
        delete _data[_data.size()-1];
        _data.pop_back();
    }
}

template <class T,int IS,class R>
void  HArray1D<T,IS,R>::erase(iterator a,iterator b)
{

}

template <class T,int IS,class R>
void  HArray1D<T,IS,R>::erase(const_iterator a,const_iterator b)
{

}

template <class T,int IS,class R>
void HArray1D<T,IS,R>::reset(const T& val)
{
    UInt sz=_data.size();
    for(i=0;i<sz;i++)
    {
        R::Memory_init(_data[i]->data,IS,val);
    }
}

template <class T,int IS,class R>
void  HArray1D<T,IS,R>::reserve(UInt sz,const T& val)
{
    int i,num=sz/IS+1;
    if(sz>0 && sz%IS==0)
        num-=1;
    UInt physize=num*IS;
    if(physize==_physize) return;
    else if(_physize>physize)
    {
        int i0=physize/IS,num=_physize/IS;
        for(i=i0;i<num;i++)
        {
            delete _data[i];
        }
        _data.resize(i0);
    }else
    {
        int num=physize/IS,i0=_physize/IS;
        for(i=i0;i<num;i++)
        {
            ATuple* anew=SL_NEW ATuple();
            R::Memory_init(anew->data,IS,val);
            _data.push_back(anew);
        }
    }
    _physize=physize;
}

template <class T,int IS,class R>
void HArray1D<T,IS,R>::resize(UInt sz,const T& val)
{
    reserve(sz,val);
    _size=sz;
}


SL_MSH_NS_END
#endif
