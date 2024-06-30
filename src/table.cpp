#include "table.hpp"

void DynTable::print_info()
{
  // update_size_hist();
  std::cout << "size: " << nkmers << "\t";
  // for (auto kv : size_hist) {
  //   std::cout << "H(" << kv.first << ")=" << kv.second << "/";
  // }
  std::cout << std::endl;
}

void DynTable::clear_rows()
{
  table.clear();
  size_hist.clear();
  nkmers = 0;
}

void DynTable::sort_columns()
{
  for (uint32_t i = 0; i < table.size(); ++i) {
    if (!table[i].empty()) {
      std::sort(table[i].begin(), table[i].end(), comp_encoding);
    }
  }
}

void DynTable::ensure_sorted_columns()
{
  for (uint32_t i = 0; i < table.size(); ++i) {
    if ((!table[i].empty()) && !std::is_sorted(table[i].begin(), table[i].end(), comp_encoding)) {
      std::sort(table[i].begin(), table[i].end(), comp_encoding);
    }
  }
}

void DynTable::update_nkmers()
{
  nkmers = 0;
  for (uint32_t i = 0; i < table.size(); ++i) {
    nkmers += table[i].size();
  }
}

void DynTable::update_size_hist()
{
  size_hist.clear();
  nkmers = 0;
  for (uint32_t i = 0; i < table.size(); ++i) {
    size_hist[table[i].size()]++;
    nkmers += table[i].size();
  }
}

void DynTable::make_unique()
{
  nkmers = 0;
  for (uint32_t i = 0; i < table.size(); ++i) {
    if (!table[i].empty()) {
      table[i].erase(std::unique(table[i].begin(), table[i].end(), eq_encoding), table[i].end());
    }
    nkmers += table[i].size();
  }
}

void DynTable::prune_columns(size_t max_size)
{
  nkmers = 0;
  for (uint32_t i = 0; i < table.size(); ++i) {
    if (table[i].size() > max_size) {
      vec<mer_t> tmp_v;
      tmp_v.reserve(max_size);
      std::sample(table[i].begin(), table[i].end(), std::back_inserter(tmp_v), max_size, gen);
      table[i] = std::move(tmp_v);
    }
    nkmers += max_size;
  }
}

void DynTable::union_table(DynTable& source)
{
  assertm(nrows == source.nrows, "Two tables differ in size.");
  if (source.table.empty()) {
    return;
  } else if (table.empty()) {
    table = std::move(source.table);
    size_hist = std::move(source.size_hist);
    nkmers = source.nkmers;
    return;
  } else {
    nkmers = 0;
    for (uint32_t i = 0; i < table.size(); ++i) {
      if (!source.table[i].empty() && !table[i].empty()) {
#ifdef NONSTD_UNION
        DynTable::union_row(table[i], source.table[i], record);
#else
        DynTable::union_row(table[i], source.table[i], record, true);
#endif
      } else if (!source.table[i].empty()) {
        table[i] = std::move(source.table[i]);
      } else {
      }
      nkmers += table[i].size();
    }
  }
}

#ifndef NONSTD_UNION
void DynTable::union_row(vec<mer_t>& dest_v,
                         vec<mer_t>& source_v,
                         record_sptr_t record,
                         bool in_place)
{
  if (in_place) {
    // Merging in-place (alternative).
    dest_v.insert(dest_v.end(),
                  std::make_move_iterator(source_v.begin()),
                  std::make_move_iterator(source_v.end()));
    std::inplace_merge(dest_v.begin(),
                       std::next(dest_v.begin(), dest_v.size() - source_v.size()),
                       dest_v.end(),
                       comp_encoding);
  } else {
    // Merging with allocation.
    vec<mer_t> tmp_v;
    tmp_v.reserve(source_v.size() + dest_v.size());
    std::merge(source_v.begin(),
               source_v.end(),
               dest_v.begin(),
               dest_v.end(),
               std::back_inserter(tmp_v),
               comp_encoding);
    // tmp_v.shrink_to_fit();
    dest_v = std::move(tmp_v);
  }
  // Deal with shared k-mers.
  auto iter = dest_v.begin(), result = dest_v.begin();
  while (++iter != dest_v.end()) {
    if (result->encoding == iter->encoding) {
      // TODO: check collisions and resolve.
      // TODO: check if subset is a node and skip.
      subset_sptr_t new_subset = std::make_shared<Subset>(iter->shash, result->shash, record);
      record->add_subset(new_subset);
      result->shash += iter->shash;
    } else if (++result != iter) {
      *result = std::move(*iter);
    } else {
      result = iter;
    }
  }
  dest_v.erase(++result, dest_v.end());
}
#else
void DynTable::union_row(vec<mer_t>& dest_v, vec<mer_t>& source_v, record_sptr_t record)
{
  vec<mer_t> temp_v;
  temp_v.reserve(source_v.size() + dest_v.size());
  auto iter_d = dest_v.begin(), iter_s = source_v.begin();
  for (; iter_s != source_v.end() && iter_d != dest_v.end(); ++iter_s) {
    while (iter_d != dest_v.end() && iter_s->encoding > iter_d->encoding) {
      temp_v.push_back(std::move(*iter_d));
      iter_d++;
    }
    while (iter_d != dest_v.end() && iter_s->encoding == iter_d->encoding) {
      // TODO: check collisions and resolve.
      // TODO: check if subset is a node and skip.
      subset_sptr_t new_subset = std::make_shared<Subset>(iter_d->shash, iter_s->shash, record);
      record->add_subset(new_subset);
      iter_s->shash += iter_d->shash;
      iter_d++;
    }
    temp_v.push_back(std::move(*iter_s));
  }
  temp_v.insert(temp_v.end(), iter_d, dest_v.end());
  // temp_v.shrink_to_fit();
  dest_v = std::move(temp_v);
}
#endif

void DynTable::fill_table(refseq_sptr_t rs)
{
  table.resize(nrows);
  while (rs->read_next_seq() && rs->set_curr_seq()) {
    rs->extract_mers(table);
  }
  // update_nkmers();
  sort_columns();
  make_unique();
}
