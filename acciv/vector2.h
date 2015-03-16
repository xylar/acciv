/****************************************
 * 2D Vector Classes
 * By Will Perone (will.perone@gmail.com)
 * Original: 9-16-2002  
 * Revised: 19-11-2003
 *          18-12-2003
 *          06-06-2004
 *
 * © 2003, This code is provided "as is" and you can use it freely as long as 
 * credit is given to Will Perone in the application it is used in
 ****************************************/

// Modified by Xylar Asay-Davis 2013-08-23

#pragma once


#include <math.h>



template <typename T>
struct vector2
{
	T x, y;


	//! trivial ctor
	vector2<T>(): x(0), y(0) {}

	//! setting ctor
	vector2<T>(const T x0, const T y0): x(x0), y(y0) {}

	//! array indexing
	T &operator [](unsigned int i) 
	{   return *(&x+i);   }

	//! array indexing
	const T &operator [](unsigned int i) const 
	{	return *(&x+i);   }  

	//! function call operator
	void operator ()(const T x0, const T y0) 
	{	x= x0; y= y0;	}

	//! test for equality
	bool operator==(const vector2<T> &v)
	{	return (x==v.x && y==v.y);	}

	//! test for inequality
	bool operator!=(const vector2<T> &v)
	{	return (x!=v.x || y!=v.y);	}

	//! set to value
	const vector2<T> &operator =(const vector2<T> &v)
	{	
		x= v.x; y= v.y;			
		return *this;
	}

	//! negation
	const vector2<T> operator -(void) const
	{	return vector2<T>(-x, -y);	}

	//! addition
	const vector2<T> operator +(const vector2<T> &v) const
	{	return vector2<T>(x+v.x, y+v.y);	}

	//! subtraction
	const vector2<T> operator -(const vector2<T> &v) const
	{   return vector2<T>(x-v.x, y-v.y);	}

	//! uniform scaling
	const vector2<T> operator *(const T num) const
	{
		vector2<T> temp(*this);			
		return temp*=num;
	}

	//! uniform scaling
	const vector2<T> operator /(const T num) const
	{
		vector2<T> temp(*this);			
		return temp/=num;
	}		

	//! addition
	const vector2<T> &operator +=(const vector2<T> &v)
	{
		x+=v.x;	y+=v.y;						
		return *this;
	}

	//! subtraction
	const vector2<T> &operator -=(const vector2<T> &v) 
	{ 
		x-=v.x;	y-=v.y;						
		return *this;	
	}

	//! uniform scaling
	const vector2<T> &operator *=(const T num)
	{
		x*=num;	y*=num;										
		return *this;
	}

	//! uniform scaling
	const vector2<T> &operator /=(const T num)
	{
		x/=num;	y/=num;										
		return *this;
	}

	//! dot product
	T dot(const vector2<T> &v) const
	{	return x*v.x + y*v.y;	}	

	//! I reluctantly overload the product of two vectors to be the dot product
	//! as this is useful for computing the L2 norm without knowing the data type
	T operator *(const vector2<T> &v) const
	{
		return this->dot(v);
	}
};

template<class T>
u_inline vector2<T> operator*(const T num,const vector2<T> & v)
{
	return vector2<T>(num*v.x,num*v.y);
}



// macro to make creating the ctors for derived vectors easier
#define VECTOR2_CTORS(name, type)   \
	/* trivial ctor */				\
	name(): vector2<type>((type)0, (type)0) {}						\
	/* down casting ctor */			\
	name(const vector2<type> &v): vector2<type>(v.x, v.y) {}	\
	/* down casting scalar ctor */			\
	name(type scalar): vector2<type>(scalar, scalar) {}	\
	/* make x,y combination into a vector */					\
	name(type x0, type y0): vector2<type>(x0, y0) {}



struct vector2i: public vector2<int>
{
	VECTOR2_CTORS(vector2i, int)
};


struct vector2ui: public vector2<unsigned int>
{
	VECTOR2_CTORS(vector2ui, unsigned int)
};


struct vector2f: public vector2<float>
{
	VECTOR2_CTORS(vector2f, float)

	//! gets the length of this vector squared
	float length_squared() const
	{	return (float)(this->dot(*this));   }

	//! gets the length of this vector
	float length() const
	{	return (float)sqrt(length_squared());   }

	//! normalizes this vector
	void normalize()
	{	*this/=length();	}

	//! returns the normalized vector
	vector2f normalized() const
	{   return  *this/length();  }

	//! reflects this vector about n
	void reflect(const vector2f &n)
	{
		vector2f orig(*this);
		project(n);
		*this= *this*2 - orig;
	}

	//! projects this vector onto v
	void project(const vector2f &v)
	{	*this=projected(v);	}

	//! returns this vector projected onto v
	vector2f projected(const vector2f &v)
	{   return v * (this->dot(v))/(v.dot(v));	}

	//! computes the angle between 2 arbitrary vectors
	static inline float angle(const vector2f &v1, const vector2f &v2)
	{   return acosf((v1.dot(v2)) / (v1.length()*v2.length()));  }

	//! computes the angle between 2 normalized arbitrary vectors
	static inline float angle_normalized(const vector2f &v1, const vector2f &v2)
	{   return acosf(v1.dot(v2));  }
};

// Added by Xylar Asay-Davis 2013-08-23
struct vector2d: public vector2<double>
{
	VECTOR2_CTORS(vector2d, double)

	//! gets the length of this vector squared
	double length_squared() const
	{	return (double)(this->dot(*this));   }

	//! gets the length of this vector
	double length() const
	{	return (double)sqrt(length_squared());   }

	//! normalizes this vector
	void normalize()
	{	*this/=length();	}

	//! returns the normalized vector
	vector2d normalized() const
	{   return  *this/length();  }

	//! reflects this vector about n
	void reflect(const vector2d &n)
	{
		vector2d orig(*this);
		project(n);
		*this= *this*2 - orig;
	}

	//! projects this vector onto v
	void project(const vector2d &v)
	{	*this=projected(v);	}

	//! returns this vector projected onto v
	vector2d projected(const vector2d &v)
	{   return v * (this->dot(v))/(v.dot(v));	}

	//! computes the angle between 2 arbitrary vectors
	static inline double angle(const vector2d &v1, const vector2d &v2)
	{   return acos((v1.dot(v2)) / (v1.length()*v2.length()));  }

	//! computes the angle between 2 normalized arbitrary vectors
	static inline double angle_normalized(const vector2d &v1, const vector2d &v2)
	{   return acos(v1.dot(v2));  }
};
