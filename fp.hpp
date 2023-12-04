#pragma once
#ifndef __FP_HPP__
#define __FP_HPP__
#include<array>

template <class T>
struct vec6
{
	T xx, yy, zz, xy, xz, yz;

	vec6(T xx, T yy, T zz, T xy, T xz, T yz)
		:
		xx{ xx },
		yy{ yy },
		zz{ zz },
		xy{ xy },
		xz{ xz },
		yz{ yz }
	{}

	vec6(T a)
		:
		vec6(T{ a }, T{ a }, T{ a }, T{ a }, T{ a }, T{ a })
	{}

	vec6() = default;

	vec6& operator+=(vec6 const& v)
	{
		xx += v.xx;
		yy += v.yy;
		zz += v.zz;
		xy += v.xy;
		xz += v.xz;
		yz += v.yz;

		return *this;
	}

	vec6 operator*(T a) const {
		return vec6(a * this->xx, a * this->yy, a * this->zz, a * this->xy, a * this->xz, a * this->yz);
	}

	vec6 operator/(T a) const {
		return vec6(this->xx / a, this->yy / a, this->zz / a, this->xy / a, this->xz / a, this->yz / a);
	}

	vec6 operator+(vec6 const& v) const {
		return vec6(this->xx + v.xx, this->yy + v.yy, this->zz + v.zz, this->xy + v.xy, this->xz + v.xz, this->yz + v.yz);
	}

	vec6 operator-(vec6 const& v) const {
		return vec6(this->xx - v.xx, this->yy - v.yy, this->zz - v.zz, this->xy - v.xy, this->xz - v.xz, this->yz - v.yz);
	}
};


template <class T>
struct vec3
{
	T x, y, z;

	vec3(T x, T y, T z) :
		x{ x },
		y{ y },
		z{ z }
	{}

	vec3(T a) :
		vec3(T{ a }, T{ a }, T{ a })
	{}

	vec3() = default;


	vec3& operator+=(vec3 const& v) {
		x += v.x;
		y += v.y;
		z += v.z;

		return *this;
	}

	vec3& operator-=(vec3 const& v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;

		return *this;
	}

	T operator*(vec3 const& v) const {
		return v.x * this->x + v.y * this->y + v.z * this->z;
	}

	vec3 operator*(T a) const {
		return vec3(a * this->x, a * this->y, a * this->z);
	}

	vec3 operator/(T a) const {
		return vec3(this->x / a, this->y / a, this->z / a);
	}

	vec3 operator+(vec3 const& v) const {
		return vec3(this->x + v.x, this->y + v.y, this->z + v.z);
	}

	vec3 operator-(vec3 const& v) const {
		return vec3(this->x - v.x, this->y - v.y, this->z - v.z);
	}
};


template<class T>
vec3<T> cross(vec3<T> const& v1, vec3<T> const& v2) {
	return vec3<T>(v1.z * v2.x - v1.x * v2.z, v1.y * v2.z - v1.z * v2.y, v1.x * v2.y - v1.y * v2.x);
}


namespace fp
{
	using fp_t = float;
	using fp3_t = vec3<fp_t>;
	using fp6_t = vec6<fp_t>;
	using fp3x3_t = std::array<std::array<fp_t, 3>, 3>;
	using fp3x3x3x3_t = std::array<std::array<std::array<std::array<fp_t, 3>, 3>, 3>, 3>;
};

#endif // !__FP_HPP__