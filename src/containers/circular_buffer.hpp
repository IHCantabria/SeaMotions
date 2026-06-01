/*
 * Copyright (c) 2025 Sergio Fernández Ruano / IHCantabria
 *
 * This file is part of SeaMotionsTimeDev.
 *
 * SeaMotions is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SeaMotions is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <cassert>
#include <utility>
#include <vector>

/**
 * @brief Fixed-capacity circular buffer with "newest-first" indexing.
 *
 * Elements are pushed via push_front(), which inserts the new item as the
 * most-recent entry (logical index 0).  Older entries shift to higher logical
 * indices automatically.  When the buffer is full the oldest entry is silently
 * overwritten, so the buffer never grows beyond the reserved capacity.
 *
 * This avoids the O(N) memory copies that would be required by
 * std::vector::insert(begin(), ...) on every time step.
 *
 * Typical usage for a radiation-memory (Duhamel) history:
 * @code
 *   CircularBuffer<std::vector<cusfloat>> hist;
 *   hist.reserve( max_steps );
 *
 *   // Each time step:
 *   hist.push_front( sigma_current );       // O(1) — no copy of existing data
 *
 *   // Duhamel loop:
 *   for ( int k = 0; k < hist.size(); ++k )
 *       ... use hist[k] ...                 // hist[0] = most recent
 * @endcode
 *
 * @tparam T  Element type.  Must be default-constructible and movable.
 */
template<typename T>
class CircularBuffer
{
public:
    CircularBuffer( ) = default;

    /**
     * @brief Pre-allocate storage for @p capacity elements.
     *
     * Safe to call more than once; subsequent calls resize the internal
     * storage and reset the buffer to an empty state.
     *
     * @param capacity  Maximum number of elements to retain.
     */
    void reserve( int capacity )
    {
        assert( capacity > 0 );
        _data.assign( capacity, T{} );
        _cap  = capacity;
        _head = 0;
        _size = 0;
    }

    /**
     * @brief Insert @p item as the most-recent element (logical index 0).
     *
     * If the buffer is already full the oldest element is overwritten.
     * Complexity: O(1) — only the internal head pointer is updated.
     *
     * @param item  New element to store (moved into the buffer).
     */
    void push_front( T item )
    {
        // Move head backward (wraps around) and write there.
        _head = ( _head - 1 + _cap ) % _cap;
        _data[_head] = std::move( item );
        if ( _size < _cap ) { ++_size; }
    }

    /**
     * @brief Access element at logical index @p k (read-only).
     *
     * Index 0 returns the most recently pushed element; index 1 returns the
     * one before that, and so on up to size()-1.
     *
     * @param k  Logical index in [0, size()).
     */
    const T& operator[]( int k ) const
    {
        assert( k >= 0 && k < _size );
        return _data[ ( _head + k ) % _cap ];
    }

    /**
     * @brief Access element at logical index @p k (mutable).
     */
    T& operator[]( int k )
    {
        assert( k >= 0 && k < _size );
        return _data[ ( _head + k ) % _cap ];
    }

    /// Number of valid elements currently stored (at most capacity()).
    int size( )     const { return _size; }

    /// Maximum number of elements the buffer can hold.
    int capacity( ) const { return _cap;  }

    /// Returns true when no elements have been pushed yet.
    bool empty( )   const { return _size == 0; }

private:
    std::vector<T>  _data;          ///< Flat storage ring
    int             _head = 0;      ///< Physical index of the most-recent element
    int             _size = 0;      ///< Number of valid elements (<= _cap)
    int             _cap  = 0;      ///< Reserved capacity
};
