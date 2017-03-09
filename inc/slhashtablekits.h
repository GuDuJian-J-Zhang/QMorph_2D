//SLHashTablekits.h
#ifndef _TPTHASHTABLEKITS_H_
#define _TPTHASHTABLEKITS_H_
#include <utility>
#include <algorithm>
#include <functional>

template <class Key> struct hash
{
};

inline size_t __stl_hash_string(const char* s)
{
    unsigned long h = 0; 
    for ( ; *s; ++s )
        h = 5*h + *s;

    return size_t(h);
}

template<>
struct hash<char*>
{
    size_t operator()(const char* s) const
    {
        return __stl_hash_string(s);
    }
};

template<>
struct hash<const char*>
{
    size_t operator()(const char* s) const
    {
        return __stl_hash_string(s);
    }
};

template<>
struct hash<char>
{
    size_t operator()(char x) const
    {
        return x;
    }
};

template<>
struct hash<unsigned char>
{
    size_t operator()(unsigned char x) const
    {
        return x;
    }
};

template<>
struct hash<signed char>
{
    size_t operator()(unsigned char x) const
    {
        return x;
    }
};

template<>
struct hash<short>
{
    size_t operator()(short x) const
    {
        return x;
    }
};

template<>
struct hash<unsigned short>
{
    size_t operator()(unsigned short x) const
    {
        return x;
    }
};

template<>
struct hash<int>
{
    size_t operator()(int x) const
    {
        return x;
    }
};

template<>
struct hash<long>
{
    size_t operator()(long x) const
    {
        return x;
    }
};

template<>
struct hash<const void*>
{
    size_t operator()(const void* ptr) const
    {
        return size_t(ptr);
    }
};

template<>
struct hash<void*>
{
    size_t operator()(void* ptr) const
    {
        return size_t(ptr);
    }
};

template<>
struct hash<size_t>
{
    size_t operator()(size_t ptr) const
    {
        return size_t(ptr);
    }
};

template <class Val>
struct hashtable_node
{
    hashtable_node()
        :_next(0)
    {
    }

    hashtable_node(const Val& val)
        :_next(0),
        _val(val)
    {
    }

    hashtable_node* _next;
    Val _val;
}; 


template <class T>
struct identity : public std::unary_function<T, T> {
    const T& operator()(const T& x) const { return x; }
};

// Note: assumes long is at least 32 bits.
enum { __stl_num_primes = 28 };

static const unsigned long __stl_prime_list[__stl_num_primes] =
{
    53ul,         97ul,         193ul,       389ul,       769ul,
    1543ul,       3079ul,       6151ul,      12289ul,     24593ul,
    49157ul,      98317ul,      196613ul,    393241ul,    786433ul,
    1572869ul,    3145739ul,    6291469ul,   12582917ul,  25165843ul,
    50331653ul,   100663319ul,  201326611ul, 402653189ul, 805306457ul, 
    1610612741ul, 3221225473ul, 4294967291ul
};

inline unsigned long NextPrime(unsigned long n)
{
    const unsigned long* first = __stl_prime_list;
    const unsigned long* last = __stl_prime_list + (int)__stl_num_primes;
    const unsigned long* pos = std::lower_bound(first, last, n);
    return pos == last ? *(last - 1) : *pos;
}


template <class Val, class Key, class HashFcn,
class ExtractKey, class EqualKey>
class SLHashTable;

template <class Val, class Key, class HashFcn,class ExtractKey, class EqualKey>
class SLHashTable
{
public:
    typedef hashtable_node<Val> HashtableNode;
    typedef size_t            size_type;
    typedef Val               value_type;
    typedef Key               key_type;
    typedef HashFcn           hasher;
    typedef EqualKey          key_equal;
    typedef ExtractKey        extract_key;
    typedef int               difference_type;
    typedef value_type*       pointer;
    typedef const value_type* const_pointer;
    typedef const value_type& const_reference;
    typedef value_type& reference;
    struct iterator 
    {
        HashtableNode* _cur;
        SLHashTable* _ht;

        iterator(HashtableNode* n, SLHashTable* tab) 
            : _cur(n), _ht(tab) {}

        iterator() {}
        reference operator*() const { return _cur->_val; }
        pointer operator->() const { return &(operator*()); }
        iterator& operator++();

        iterator operator++(int)
        {
            iterator tmp = *this;
            ++*this;
            return tmp;
        }

        bool operator==(const iterator& it) const
        { return _cur == it._cur; }
        bool operator!=(const iterator& it) const
        { return _cur != it._cur; }
    };


    struct const_iterator 
    {
        const HashtableNode* _cur;
        const SLHashTable* _ht;

        const_iterator(const HashtableNode* n, const SLHashTable* tab)
            : _cur(n), _ht(tab) {}

        const_iterator() {}
        const_iterator(const iterator& it) 
            : _cur(it._cur), _ht(it._ht) {}
        const_reference operator*() const { return _cur->_val; }
        const_pointer operator->() const { return &(operator*()); }
        const_iterator& operator++();

        const_iterator operator++(int)
        {
            const_iterator tmp = *this;
            ++*this;
            return tmp;
        }

        bool operator==(const const_iterator& it) const 
        { return _cur == it._cur; }
        bool operator!=(const const_iterator& it) const 
        { return _cur != it._cur; }
    };

    struct allocator_type{};
    friend struct const_iterator;
    friend struct iterator;
protected:
    std::vector<HashtableNode*>   _buckets;
    size_type                     _num;
    hasher                        _hash;
    key_equal                     _equals;
    extract_key                   _get_key;
public:
    SLHashTable(size_type n,
        const hasher&    hf,
        const key_equal&   eql,
        const extract_key& ext)
        :_hash(hf),
        _equals(eql),
        _get_key(ext),
        _num(0)
    {
        initialize_buckets(n);
    }

    SLHashTable(size_type n,
        const hasher&    hf,
        const key_equal&   eql)
        :_hash(hf),
        _equals(eql),
        _get_key(extract_key()),
        _num(0)
    {
        initialize_buckets(n);
    }

    ~SLHashTable() { clear(); }

    SLHashTable& operator= (const SLHashTable& ht)
    {
        if (&ht != this) 
        {
            clear();
            _hash = ht._hash;
            _equals = ht._equals;
            _get_key = ht._get_key;
            copy_from(ht);
        }
        return *this;
    }

    void swap(SLHashTable& ht)
    {
        std::swap(_hash, ht._hash);
        std::swap(_equals, ht._equals);
        std::swap(_get_key, ht._get_key);
        _buckets.swap(ht._buckets);
        std::swap(_num, ht._num);
    }

    size_type size() const { return _num; }
    size_type max_size() const { return size_type(-1); }
    bool empty() const { return size() == 0; }

    iterator begin()
    { 
        for (size_type n = 0; n < _buckets.size(); ++n)
            if (_buckets[n])
                return iterator(_buckets[n], this);
        return end();
    }

    iterator end() { return iterator(0, this); }

    const_iterator begin() const
    {
        for (size_type n = 0; n < _buckets.size(); ++n)
            if (_buckets[n])
                return const_iterator(_buckets[n], this);
        return end();
    }

    const_iterator end() const { return const_iterator(0, this); }

    std::pair<iterator, bool> insert_unique(const value_type& obj)
    {
        resize(_num + 1);
        return insert_unique_noresize(obj);
    }

    size_type erase(const key_type& key)
    {
        const size_type n = bkt_num_key(key);
        HashtableNode* first = _buckets[n];
        size_type erased = 0;

        if (first) 
        {
            HashtableNode* cur = first;
            HashtableNode* next = cur->_next;
            while (next) 
            {
                if (_equals(_get_key(next->_val), key)) 
                {
                    cur->_next = next->_next;
                    delete_node(next);
                    next = cur->_next;
                    ++erased;
                    --_num;
                }else 
                {
                    cur = next;
                    next = cur->_next;
                }
            }
            if (_equals(_get_key(first->_val), key)) 
            {
                _buckets[n] = first->_next;
                delete_node(first);
                ++erased;
                --_num;
            }
        }
        return erased;
    }

    void erase(const iterator& it)
    {
        HashtableNode* p = it._cur;
        if (p) 
        {
            const size_type n = bkt_num(p->_val);
            HashtableNode* cur = _buckets[n];

            if (cur == p) 
            {
                _buckets[n] = cur->_next;
                delete_node(cur);
                --_num;
            }else 
            {
                HashtableNode* next = cur->_next;
                while (next) 
                {
                    if (next == p) 
                    {
                        cur->_next = next->_next;
                        delete_node(next);
                        --_num;
                        break;
                    } else 
                    {
                        cur = next;
                        next = cur->_next;
                    }
                }
            }
        }
    }

    void erase(const const_iterator& it)
    {
        erase(iterator(const_cast<HashtableNode*>(it._cur),const_cast<SLHashTable*>(it._ht)));
    }

    iterator find(const key_type& key) 
    {
        size_type n = bkt_num_key(key);
        HashtableNode* first;
        for ( first = _buckets[n];
            first && !_equals(_get_key(first->_val), key);
            first = first->_next)
        {}
        return iterator(first, this);
    } 

    const_iterator find(const key_type& key) const
    {
        size_type n = bkt_num_key(key);
        const HashtableNode* first;
        for ( first = _buckets[n];
            first && !_equals(_get_key(first->_val), key);
            first = first->_next)
        {}
        return const_iterator(first, this);
    } 

    size_type count(const key_type& key) const
    {
        const size_type n = bkt_num_key(key);
        size_type __result = 0;

        for (const HashtableNode* cur = _buckets[n]; cur; cur = cur->_next)
            if (_equals(_get_key(cur->_val), key))
                ++__result;
        return __result;
    }

    void insert_unique(const value_type* f, const value_type* l)
    {
        size_type n = l - f;
        resize(_num + n);
        for ( ; n > 0; --n, ++f)
            insert_unique_noresize(*f);
    }

    void insert_equal(const value_type* f, const value_type* l)
    {
        size_type n = l - f;
        resize(_num + n);
        for ( ; n > 0; --n, ++f)
            insert_equal_noresize(*f);
    }

    void insert_unique(const_iterator f, const_iterator l)
    {
        size_type n = 0;
        distance(f, l, n);
        resize(_num + n);
        for ( ; n > 0; --n, ++f)
            insert_unique_noresize(*f);
    }

    void insert_equal(const_iterator f, const_iterator l)
    {
        size_type n = 0;
        distance(f, l, n);
        resize(_num + n);
        for ( ; n > 0; --n, ++f)
            insert_equal_noresize(*f);
    }

    size_type bucket_count() const { return _buckets.size(); }

    size_type max_bucket_count() const
    { return __stl_prime_list[(int)__stl_num_primes - 1]; } 

    size_type elems_in_bucket(size_type __bucket) const
    {
        size_type __result = 0;
        for (HashtableNode* cur = _buckets[__bucket]; cur; cur = cur->_next)
            __result += 1;
        return __result;
    }


    iterator insert_equal(const value_type& obj)
    {
        resize(_num + 1);
        return insert_equal_noresize(obj);
    }

    std::pair<iterator, bool> insert_unique_noresize(const value_type& obj)
    {
        const size_type n = bkt_num(obj);
        HashtableNode* first = _buckets[n];

        for (HashtableNode* cur = first; cur; cur = cur->_next) 
            if (_equals(_get_key(cur->_val), _get_key(obj)))
                return std::pair<iterator, bool>(iterator(cur, this), false);

        HashtableNode* tmp = new_node(obj);
        tmp->_next = first;
        _buckets[n] = tmp;
        ++_num;
        return std::pair<iterator, bool>(iterator(tmp, this), true);
    }

    iterator insert_equal_noresize(const value_type& obj)
    {
        const size_type n = bkt_num(obj);
        HashtableNode* first = _buckets[n];

        for (HashtableNode* cur = first; cur; cur = cur->_next) 
            if (_equals(_get_key(cur->_val), _get_key(obj))) {
                HashtableNode* tmp = new_node(obj);
                tmp->_next = cur->_next;
                cur->_next = tmp;
                ++_num;
                return iterator(tmp, this);
            }

            HashtableNode* tmp = new_node(obj);
            tmp->_next = first;
            _buckets[n] = tmp;
            ++_num;
            return iterator(tmp, this);
    }

    reference find_or_insert(const value_type& obj)
    {
        resize(_num + 1);

        size_type n = bkt_num(obj);
        HashtableNode* first = _buckets[n];

        for (HashtableNode* cur = first; cur; cur = cur->_next)
            if (_equals(_get_key(cur->_val), _get_key(obj)))
                return cur->_val;

        HashtableNode* tmp = new_node(obj);
        tmp->_next = first;
        _buckets[n] = tmp;
        ++_num;
        return tmp->_val;
    }

    std::pair<iterator,iterator> equal_range(const key_type& key)
    {
        typedef std::pair<iterator, iterator> pii;
        const size_type n = bkt_num_key(key);

        for (HashtableNode* first = _buckets[n]; first; first = first->_next)
            if (_equals(_get_key(first->_val), key)) 
            {
                for (HashtableNode* cur = first->_next; cur; cur = cur->_next)
                    if (!_equals(_get_key(cur->_val), key))
                        return pii(iterator(first, this), iterator(cur, this));
                for (size_type m = n + 1; m < _buckets.size(); ++m)
                    if (_buckets[m])
                        return pii(iterator(first, this),
                        iterator(_buckets[m], this));
                return pii(iterator(first, this), end());
            }
            return pii(end(), end());
    }

    std::pair<const_iterator, const_iterator> equal_range(const key_type& key) const
    {
        typedef std::pair<const_iterator, const_iterator> pii;
        const size_type n = bkt_num_key(key);

        for (const HashtableNode* first = _buckets[n] ;
            first; 
            first = first->_next) 
        {
            if (_equals(_get_key(first->_val), key)) 
            {
                for (const HashtableNode* cur = first->_next;
                    cur;
                    cur = cur->_next)
                    if (!_equals(_get_key(cur->_val), key))
                        return pii(const_iterator(first, this),
                        const_iterator(cur, this));
                for (size_type m = n + 1; m < _buckets.size(); ++m)
                    if (_buckets[m])
                        return pii(const_iterator(first, this),
                        const_iterator(_buckets[m], this));
                return pii(const_iterator(first, this), end());
            }
        }
        return pii(end(), end());
    }



    void clear()
    {
        for (size_type i = 0; i < _buckets.size(); ++i) 
        {
            HashtableNode* cur = _buckets[i];
            while (cur != 0) 
            {
                HashtableNode* next = cur->_next;
                delete_node(cur);
                cur = next;
            }
            _buckets[i] = 0;
        }
        _num = 0;
    }


protected:
    size_type next_size(size_type n) const{ return NextPrime(n); }

    size_type bkt_num_key(const key_type& key) const
    {
        return bkt_num_key(key, _buckets.size());
    }

    size_type bkt_num(const value_type& obj) const
    {
        return bkt_num_key(_get_key(obj));
    }

    size_type bkt_num_key(const key_type&  key, size_type n) const
    {
        return _hash(key) % n;
    }

    size_type bkt_num(const value_type& obj, size_type n) const
    {
        return bkt_num_key(_get_key(obj), n);
    }

    void initialize_buckets(size_type n)
    {
        const size_type bn = next_size(n);
        _buckets.reserve(bn);
        _buckets.insert(_buckets.end(), bn, (HashtableNode*) 0);
        _num = 0;
    }

    HashtableNode* new_node(const value_type& obj)
    {
        HashtableNode* n = new HashtableNode(obj);
        n->_next = 0;
        //n->_val=obj;
        return n;
    }

    void delete_node(HashtableNode* n)
    {
        delete n;
    }


    void copy_from(const SLHashTable& ht)
    {
        _buckets.clear();
        _buckets.reserve(ht._buckets.size());
        _buckets.insert(_buckets.end(), ht._buckets.size(), (HashtableNode*) 0);
        try 
        {
            for (size_type i = 0; i < ht._buckets.size(); ++i) 
            {
                const HashtableNode* cur = ht._buckets[i];
                if (cur) 
                {
                    HashtableNode* copy = new_node(cur->_val);
                    _buckets[i] = copy;
                    HashtableNode* next = 0;
                    for (next = cur->_next;next;cur = next,next = cur->_next) 
                    {
                        copy->_next = new_node(next->_val);
                        copy = copy->_next;
                    }
                }
            }
            _num = ht._num;
        }catch(...)
        {
            clear();
        }
    }


    void resize(size_type num)
    {
        const size_type old_n = _buckets.size();
        if (num > old_n) 
        {
            const size_type n = next_size(num);
            if (n > old_n) 
            {
                std::vector<HashtableNode*> tmp(n, (HashtableNode*)(0),_buckets.get_allocator());
                try 
                {
                    for (size_type bucket = 0; bucket < old_n; ++bucket) 
                    {
                        HashtableNode* first = _buckets[bucket];
                        while (first) 
                        {
                            size_type new_bucket = bkt_num(first->_val, n);
                            _buckets[bucket] = first->_next;
                            first->_next = tmp[new_bucket];
                            tmp[new_bucket] = first;
                            first = _buckets[bucket];          
                        }
                    }
                    _buckets.swap(tmp);
                }catch(...) 
                {
                    for (size_type bucket = 0; bucket < tmp.size(); ++bucket) 
                    {
                        while (tmp[bucket]) 
                        {
                            HashtableNode* next = tmp[bucket]->_next;
                            delete_node(tmp[bucket]);
                            tmp[bucket] = next;
                        }
                    }
                    throw;
                }
            }
        }
    }

    void erase_bucket(const size_type n, HashtableNode* last)
    {
        HashtableNode* cur = _buckets[n];
        while (cur != last) 
        {
            HashtableNode* next = cur->_next;
            delete_node(cur);
            cur = next;
            _buckets[n] = cur;
            --_num;
        }
    }

    void erase_bucket(const size_type n, HashtableNode* first, HashtableNode* last)
    {
        HashtableNode* cur = _buckets[n];
        if (cur == first)
            erase_bucket(n, last);
        else 
        {
            HashtableNode* next;
            for (next = cur->_next; 
                next != first; 
                cur = next, next = cur->_next)
                ;
            while (next != last) 
            {
                cur->_next = next->_next;
                delete_node(next);
                next = cur->_next;
                --_num;
            }
        }
    }

};


template <class Val, class Key, class HashFcn,class ExtractKey, class EqualKey>
typename SLHashTable<Val,Key,HashFcn,ExtractKey,EqualKey>::const_iterator&
SLHashTable<Val,Key,HashFcn,ExtractKey,EqualKey>::const_iterator::operator++()
{
    const HashtableNode* old = _cur;
    _cur = _cur->_next;
    if (!_cur) {
        size_type bucket = _ht->bkt_num(old->_val);
        while (!_cur && ++bucket < _ht->_buckets.size())
            _cur = _ht->_buckets[bucket];
    }
    return *this;
}


template <class Val, class Key, class HashFcn,class ExtractKey, class EqualKey>
typename SLHashTable<Val,Key,HashFcn,ExtractKey,EqualKey>::iterator&
SLHashTable<Val,Key,HashFcn,ExtractKey,EqualKey>::iterator::operator++()
{
    const HashtableNode* old = _cur;
    _cur = _cur->_next;
    if (!_cur) 
    {
        size_type bucket = _ht->bkt_num(old->_val);
        while (!_cur && ++bucket < _ht->_buckets.size())
            _cur = _ht->_buckets[bucket];
    }
    return *this;
}

#endif
