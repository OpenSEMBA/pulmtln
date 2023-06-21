#pragma once

#include <cstddef>
#include <iostream>
#include <string>
#include <sstream>

namespace semba {
namespace util {

template <typename T>
class Identification {
    template <typename I>
    friend std::istream& operator>>(std::istream&  input,
                                    Identification<I>& id);
    template <typename I>
    friend std::ostream& operator<<(std::ostream& output,
                                    const Identification<I>& id);
public:
    Identification() = default;
    explicit Identification(const std::size_t id);
    Identification(const Identification& rhs);
    virtual ~Identification() = default;

    Identification& operator =(const Identification& rhs);
    Identification& operator+=(const Identification& rhs);

    bool operator==(const Identification& rhs) const;
    bool operator!=(const Identification& rhs) const;
    bool operator< (const Identification& rhs) const;
    bool operator<=(const Identification& rhs) const;
    bool operator> (const Identification& rhs) const;
    bool operator>=(const Identification& rhs) const;

    Identification  operator+ (const Identification& rhs) const;
    Identification& operator++();
    Identification  operator++(int);

    std::size_t toInt() const;
    std::string toStr() const;

private:
    std::size_t id_{ 0 };
};


template<typename T>
Identification<T>::Identification(const std::size_t id)
    : id_(id) {

}

template<typename T>
Identification<T>::Identification(const Identification& rhs)
    : id_(rhs.id_) {

}

template<typename T>
Identification<T>& Identification<T>::operator=(const Identification& rhs) {
    if (this == &rhs) {
        return *this;
    }
    id_ = rhs.id_;
    return *this;
}

template<typename T>
Identification<T>& Identification<T>::operator+=(const Identification& rhs) {
    id_ += rhs.id_;
    return *this;
}

template<typename T>
bool Identification<T>::operator==(const Identification& rhs) const {
    return id_ == rhs.id_;
}

template<typename T>
bool Identification<T>::operator!=(const Identification& rhs) const {
    return id_ != rhs.id_;
}

template<typename T>
bool Identification<T>::operator<(const Identification& rhs) const {
    return id_ < rhs.id_;
}

template<typename T>
bool Identification<T>::operator<=(const Identification& rhs) const {
    return id_ <= rhs.id_;
}

template<typename T>
bool Identification<T>::operator>(const Identification& rhs) const {
    return id_ > rhs.id_;
}

template<typename T>
bool Identification<T>::operator>=(const Identification& rhs) const {
    return id_ >= rhs.id_;
}

template<typename T>
Identification<T> Identification<T>::operator+(
    const Identification& rhs) const {
    return Identification(id_ + rhs.id_);
}

template<typename T>
Identification<T>& Identification<T>::operator++() {
    id_++;
    return *this;
}

template<typename T>
Identification<T> Identification<T>::operator++(int) {
    Identification copy(*this);
    id_++;
    return copy;
}

template<typename T>
std::size_t Identification<T>::toInt() const {
    return id_;
}

template<typename T>
std::string Identification<T>::toStr() const {
    std::stringstream aux;
    aux << id_;
    return aux.str();
}

template <typename I>
std::istream& operator>>(std::istream& input, Identification<I>& id) {
    input >> id.id_;
    return input;
}

template <typename I>
std::ostream& operator<<(std::ostream& output, const Identification<I>& id) {
    output << id.id_;
    return output;
}

}
} 