#include "table.hpp"

#define assertm(exp, msg) assert(((void)msg, exp))

void DynTable::clear_rows() {
  table.clear();
  table.resize(nrows);
  nkmers = 0;
  size_hist.clear();
}

void DynTable::sort_columns() {
  for (uint32_t i = 0; i < nrows; ++i) {
    if (!table[i].empty()) {
      std::sort(table[i].begin(), table[i].end(), comp_encoding);
    }
  }
}

void DynTable::ensure_sorted_columns() {
  for (uint32_t i = 0; i < nrows; ++i) {
    if ((!table[i].empty()) &&
        !std::is_sorted(table[i].begin(), table[i].end(), comp_encoding)) {
      std::sort(table[i].begin(), table[i].end(), comp_encoding);
    }
  }
}

void DynTable::update_size_hist() {
  size_hist.clear();
  nkmers = 0;
  for (uint32_t i = 0; i < nrows; ++i) {
    size_hist[table[i].size()]++;
    nkmers += table[i].size();
  }
}

void DynTable::make_unique() {
  for (uint32_t i = 0; i < nrows; ++i) {
    if (!table[i].empty()) {
      nkmers -= table[i].size();
      table[i].erase(std::unique(table[i].begin(), table[i].end(), eq_encoding),
                     table[i].end());
      nkmers += table[i].size();
    }
  }
}

void DynTable::prune_columns(size_t max_size) {
  for (uint32_t i = 0; i < nrows; ++i) {
    if (table[i].size() > max_size) {
      nkmers -= table[i].size();
      vec<mer_t> tmp_v;
      tmp_v.reserve(max_size);
      std::sample(table[i].begin(), table[i].end(), std::back_inserter(tmp_v),
                  max_size, gen);
      table[i] = std::move(tmp_v);
      nkmers += max_size;
    }
  }
}

void DynTable::union_table(DynTable &source) {
  assertm(nrows == source.nrows, "Two tables differ in size.");
  record->union_record(source.record);
  for (uint32_t i = 0; i < nrows; ++i) {
    if (!source.table[i].empty() && !table[i].empty()) {
      nkmers -= table[i].size();
#ifdef NONSTD_UNION
      DynTable::union_row(table[i], source.table[i], record);
#else
      DynTable::union_row(table[i], source.table[i], record, false);
#endif
      nkmers += table[i].size();
    } else if (!source.table[i].empty()) {
      table[i] = std::move(source.table[i]);
      nkmers += source.table[i].size();
    } else {
      continue;
    }
  }
}

#ifndef NONSTD_UNION
void DynTable::union_row(vec<mer_t> &dest_v, vec<mer_t> &source_v,
                         record_sptr_t record, bool in_place) {
  if (in_place) {
    // Merging in-place (alternative).
    dest_v.insert(dest_v.end(), std::make_move_iterator(source_v.begin()),
                  std::make_move_iterator(source_v.end()));
    std::inplace_merge(
        dest_v.begin(),
        std::next(dest_v.begin(), dest_v.size() - source_v.size()),
        dest_v.end(), comp_encoding);
  } else {
    // Merging with allocation.
    vec<mer_t> tmp_v;
    tmp_v.reserve(source_v.size() + dest_v.size());
    std::merge(source_v.begin(), source_v.end(), dest_v.begin(), dest_v.end(),
               std::back_inserter(tmp_v), comp_encoding);
    // tmp_v.shrink_to_fit();
    dest_v = std::move(tmp_v);
  }
  // Deal with shared k-mers.
  auto iter = dest_v.begin(), result = dest_v.begin();
  while (++iter != dest_v.end()) {
    if (result->encoding == iter->encoding) {
      // TODO: check collisions and resolve.
      // TODO: check if subset is a node and skip.
      subset_sptr_t new_subset =
          make_shared<Subset>(iter_d->shash, iter_s->shash, record);
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
void DynTable::union_row(vec<mer_t> &dest_v, vec<mer_t> &source_v,
                         record_sptr_t record) {
  vec<mer_t> temp_v;
  temp_v.reserve(source_v.size() + dest_v.size());
  auto iter_d = dest_v.begin(), iter_s = source_v.begin();
  for (; iter_s != source_v.end(); ++iter_s) {
    while (iter_s->encoding > iter_d->encoding) {
      temp_v.push_back(std::move(*iter_d));
      iter_d++;
    }
    while (iter_s->encoding == iter_d->encoding) {
      // TODO: check collisions and resolve.
      // TODO: check if subset is a node and skip.
      subset_sptr_t new_subset =
          make_shared<Subset>(iter_d->shash, iter_s->shash, record);
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

void DynTable::fill_table(builder_sptr_t builder) {
  while (builder->read_next_seq() && builder->set_curr_seq()) {
    builder->extract_mers(table);
  }
  sort_columns();
  make_unique();
  // TODO: check when the builder class is finalized.
}
