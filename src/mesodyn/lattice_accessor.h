#ifndef LATTICE_ACCESSOR_H
#define LATTICE_ACCESSOR_H

#include <unistd.h> //size_t
#include <map>
#include <functional>

class Lattice;

enum class Dimension {
        X,
        Y,
        Z,
        ALL
};

enum Dimensionality {
        one_D = 1,
        two_D = 2,
        three_D = 3
};

/* Do **NOT** inherit this publicly, accidental upcasting will cause a boatload of trouble */

class Lattice_accessor {
    static constexpr uint8_t SYSTEM_EDGE_OFFSET = 1;
    static constexpr uint8_t BOUNDARIES = 2;
    typedef std::map<Dimension, size_t> Coordinate;

  public:
    Lattice_accessor(const Lattice* Lat);

    const size_t MX, MY, MZ;
    const size_t system_size;
    //in lattice: gradients
    const Dimensionality dimensionality;

    const Coordinate coordinate(size_t index);

    void skip_bounds(std::function<void(size_t, size_t, size_t)> function) noexcept;

    void skip_bounds_x_major(std::function<void(size_t, size_t, size_t)> function) noexcept;

    void full_system_plus_direction_neighborlist(std::function<void(size_t, size_t, size_t)> function) noexcept;

    void full_system_minus_direction_neighborlist(std::function<void(size_t, size_t, size_t)> function) noexcept;

    void system_plus_bounds(std::function<void(size_t, size_t, size_t)> function) noexcept;

    void system_plus_bounds_x_major(std::function<void(size_t, size_t, size_t)> function) noexcept;

    size_t index (const size_t x, const size_t y, const size_t z) noexcept;

    size_t index (const size_t x, const size_t y, const size_t z) const noexcept;

    void x0_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept;

    void xm_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept;

    void y0_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept;

    void ym_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept;

    void z0_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept;

    void zm_boundary(std::function<void(size_t, size_t, size_t)> function) noexcept;

  private:
      //in lattice: JX
    const size_t jump_x;
    //in lattice: JY
    const size_t jump_y;
    //in lattice: JZ
    const size_t jump_z;
    //in lattice: M

};

#endif