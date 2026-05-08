/*
 * Copyright (c) 2025 Sergio Fernandez Ruano / IHCantabria
 *
 * This file is part of SeaMotions Software.
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

#include <algorithm>
#include <cstddef>
#include <ostream>
#include <unordered_map>
#include <string>
#include <vector>

struct MemoryReportEntry
{
    std::string name;
    std::size_t bytes = 0;
};

inline std::string memory_report_path(const std::string& prefix, const std::string& name)
{
    if (prefix.empty())
    {
        return name;
    }
    return prefix + "." + name;
}

inline void add_memory_entry(std::vector<MemoryReportEntry>& entries, const std::string& name, std::size_t bytes)
{
    if (bytes > 0)
    {
        entries.push_back({ name, bytes });
    }
}

inline std::size_t sum_memory_entries(const std::vector<MemoryReportEntry>& entries)
{
    std::size_t total = 0;
    for (const auto& entry : entries)
    {
        total += entry.bytes;
    }
    return total;
}

inline double bytes_to_mb(std::size_t bytes)
{
    return static_cast<double>( bytes ) / ( 1024.0 * 1024.0 );
}

inline std::string memory_report_root(const std::string& name)
{
    std::size_t pos = name.find( '.' );
    if ( pos == std::string::npos )
    {
        return name;
    }
    return name.substr( 0, pos );
}

inline void write_memory_report_csv(std::ostream& out, std::vector<MemoryReportEntry> entries)
{
    std::sort(
                entries.begin(),
                entries.end(),
                [](const MemoryReportEntry& a, const MemoryReportEntry& b)
                {
                    if (a.bytes != b.bytes)
                    {
                        return a.bytes > b.bytes;
                    }
                    return a.name < b.name;
                }
            );

    std::unordered_map<std::string, std::size_t> object_totals;
    for (const auto& entry : entries)
    {
        object_totals[memory_report_root(entry.name)] += entry.bytes;
    }

    out << "name,mb\n";
    for (const auto& entry : entries)
    {
        out << entry.name << "," << bytes_to_mb( entry.bytes ) << "\n";
    }
    for (const auto& item : object_totals)
    {
        out << "OBJECT_TOTAL." << item.first << "," << bytes_to_mb( item.second ) << "\n";
    }
    out << "TOTAL," << bytes_to_mb( sum_memory_entries(entries) ) << "\n";
}
