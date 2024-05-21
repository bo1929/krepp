#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

typedef uint64_t h_t;
struct esset {
  h_t ss;
  uint32_t c;
  uint32_t a;
  uint32_t ke = 0;
};

typedef std::vector<std::vector<std::pair<encT, h_t>>> vvec_merrec;
typedef std::unordered_map<std::string, h_t> um_str_h_t;
typedef std::unordered_map<h_t, std::string> um_uint32_h_t;
typedef std::unordered_map<h_t, esset> um_h_esset;

typedef struct node_t {
  std::string name;
  double blen; // branch length
  node_t *up;
  int nbranches;
  int mbranches;
  node_t **branches;
  int serial;
  h_t sh;
} node_t;

int serial = 0;

node_t *alloc() {
  node_t *nd = new node_t;
  nd->name = "";
  nd->blen = -1;
  nd->mbranches = 2;
  nd->nbranches = 0;
  nd->up = NULL;
  nd->branches = (node_t **)calloc(nd->mbranches, sizeof(node_t *));
  nd->serial = serial++;
  nd->sh = 0;
  return nd;
}

void mysplit2(char *str, std::vector<std::string> &vec) {
  int i = 0;
  char buf[strlen(str)]; // THIS IS ENOUGH
  memset(buf, 0, strlen(str));
  int at = 0;
  for (; i < strlen(str); i++) {
    if (str[i] == '\t' || str[i] == ' ' || str[i] == ';')
      continue;

    if (str[i] == '(' || str[i] == ')' || str[i] == ':' || str[i] == ',') {
      if (strlen(buf) > 0) {
        vec.push_back(buf);
        memset(buf, 0, strlen(str));
        at = 0;
      }
      vec.push_back(std::string() + str[i]);
      continue;
    } else {
      buf[at++] = str[i];
    }
  }
  if (strlen(buf) > 0) {
    vec.push_back(buf);
    memset(buf, 0, strlen(str));
    at = 0;
  }
}

int atter = 0; // dag

// production rule: prod := (prod,prod)nam:len
node_t *parse(std::vector<std::string> &vec) {

  node_t *nd = alloc();
  int seri = nd->serial;
  if (atter >= vec.size())
    return nd;
  // catch leaflike case
  if (vec[atter] != "(") {
    if (vec[atter] != ":") {
      nd->name = vec[atter];
      atter++;
    }
    if (vec[atter] == ":") {
      nd->blen = atof(vec[atter + 1].c_str());
      atter += 2;
    }
    return nd;
  }

  // catch recursive cases;
  if (vec[atter] == "(") {
    while (1) {
      atter++;
      node_t *anode = parse(vec);
      anode->up = nd;
      if (nd->nbranches == nd->mbranches) {
        nd->mbranches = 2 * nd->mbranches;
        nd->branches =
            (node_t **)realloc(nd->branches, sizeof(node_t *) * nd->mbranches);
        /* fprintf(stderr, "reallocing\n"); */
      }
      nd->branches[nd->nbranches++] = anode;
      if (vec[atter] == ",")
        continue;
      else
        break;
    }
  }
  if (vec[atter] == ")") {
    atter++;
    if (vec[atter] ==
        ")") // catch the case when there is no name or branch length
      return nd;
  }

  if (vec[atter] != ":") {
    nd->name = vec[atter];
    atter++;
  }
  if (vec[atter] == ":") {
    nd->blen = atof(vec[atter + 1].c_str());
    atter += 2;
  }
  if (!nd->name.empty() && nd->name[nd->name.length() - 1] == '\n') {
    nd->name.erase(nd->name.length() - 1);
  }
  return nd;
}

void print_node(FILE *fp, node_t *nd) {
  fprintf(fp, "----------\nnd->serial: %d this: %p\n", nd->serial, nd);
  //  fprintf(fp,"nd->left: %p nd->right: %p nd->up:
  //  %p\n",nd->left,nd->right,nd->up);
  fprintf(fp, "nd->name: %s nd->blen: %f nd->up:%p\n----------\n",
          (nd->name).c_str(), nd->blen, nd->up);
}

void serialize(node_t *nd, node_t **lst) {
  assert(nd);
  lst[nd->serial] = nd;
  for (int i = 0; i < nd->nbranches; i++)
    if (nd->branches[i])
      serialize(nd->branches[i], lst);
}

uint64_t
post_order_traversal(node_t *nd, std::vector<node_t *> &po_vec,
                     std::unordered_map<std::string, std::string> &le_fpaths) {

  h_t sh = 0;
  for (int i = 0; i < nd->nbranches; ++i) {
    sh += post_order_traversal(nd->branches[i], po_vec, le_fpaths);
  }
  if (!nd->nbranches && le_fpaths.contains(nd->name)) {
    uint32_t a1, a2;
    MurmurHash3_x86_32((nd->name).c_str(), (nd->name).length(), 0, &a1);
    MurmurHash3_x86_32((nd->name).c_str(), (nd->name).length(), 1, &a2);
    sh = sh | static_cast<uint64_t>(a1);
    sh = sh << 32;
    sh = sh | static_cast<uint64_t>(a2);
  }
  nd->sh = sh;
  po_vec.push_back(nd);
  return sh;
}

void process_genome(node_t *nd, inputHandler<encT> &pI, vvec_merrec &table,
                    um_h_esset &shash2desc) {
  if (!shash2desc.contains(nd->sh)) {
    shash2desc[nd->sh].ss = nd->sh;
    shash2desc[nd->sh].c = 0;
    shash2desc[nd->sh].a = 1;
  }
  for (auto &lsh_enc : pI.lsh_enc_vec) {
    bool seen = false;
    for (unsigned int i = 0; i < table[lsh_enc.first].size() && !seen; ++i) {
      if (lsh_enc.second == table[lsh_enc.first][i].first) {
        h_t shash = table[lsh_enc.first][i].second;
        h_t new_shash = shash + nd->sh;
        // uint64_t conc_sghash = 0;
        // conc_sghash += shash;
        // conc_sghash= (conc_sghash << 32) + nd->sh;
        // MurmurHash3_x86_32(&conc_sghash, 8, 0, &new_shash);

        shash2desc[shash].c--;
        bool existing_hash = shash2desc.contains(new_shash);

        if (existing_hash) {
          shash2desc[new_shash].c++;
        } else {
          if (shash2desc[shash].a > 1) {
            auto &m_rec = shash2desc[shash];
            h_t sh_lc = shash2desc[m_rec.ss].a < shash2desc[shash - m_rec.ss].a
                            ? m_rec.ss
                            : shash - m_rec.ss;
            h_t ss_new = sh_lc + nd->sh;
            shash2desc[new_shash].ss = ss_new;
            shash2desc[new_shash].c = 1;
            shash2desc[new_shash].a = 1 + shash2desc[shash].a;
            if (!shash2desc.contains(ss_new)) {
              shash2desc[ss_new].ss = sh_lc;
              shash2desc[ss_new].c = 0;
              shash2desc[ss_new].a = 1 + shash2desc[sh_lc].a;
            }
          } else {
            shash2desc[new_shash].ss = shash;
            shash2desc[new_shash].c = 1;
            shash2desc[new_shash].a = 1 + shash2desc[shash].a;
          }
        }
        seen = true;
        table[lsh_enc.first][i].second = new_shash;
      }
    }
    if (!seen) {
      table[lsh_enc.first].push_back(std::make_pair(lsh_enc.second, nd->sh));
      shash2desc[nd->sh].c++;
    }
  }
}
