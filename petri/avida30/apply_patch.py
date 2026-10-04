#!/usr/bin/env python3
"""Apply the material-draw edits to an Avida checkout (run from avida-core/source)."""
import re, sys

def edit(path, old, new, count=1):
    """Replace old with new; trailing whitespace on each line of old is matched loosely."""
    s = open(path).read()
    lines = old.split('\n')
    pat = '\n'.join(re.escape(l.rstrip()) + (r'[ \t]*' if i < len(lines) - 1 else '') for i, l in enumerate(lines))
    found = list(re.finditer(pat, s))
    if len(found) != count:
        sys.exit(f"{path}: expected {count} match(es), found {len(found)}: {old[:60]!r}")
    s = re.sub(pat, lambda m: new, s)
    open(path, 'w').write(s)

# ---------------- config ----------------
edit('main/cAvidaConfig.h',
'  CONFIG_ADD_GROUP(DEPRECATED_GROUP,',
'''  // -------- Conserved material (Code Petri Dish port) --------
  CONFIG_ADD_GROUP(MATERIAL_GROUP, "Conserved-material copying (h-copy-mat)");
  CONFIG_ADD_VAR(MATERIAL_MODE, int, 0, "0 = off (h-copy-mat behaves as h-copy)\\n2 = local integer material pools per cell and instruction");
  CONFIG_ADD_VAR(MATERIAL_DENSITY, double, 8.0, "Letters of material per cell at the start, uniformly random over the instruction set");
  CONFIG_ADD_VAR(MATERIAL_REACH, int, 4, "Neighbours searched after the own cell: 0, 4 (N,E,S,W) or 8");
  CONFIG_ADD_VAR(MATERIAL_DIFFUSE, double, 0.0, "Per update, cells*value attempts to move one free letter to a random orthogonal neighbour");
  CONFIG_ADD_VAR(MATERIAL_SEED, int, 0, "Seed of the material streams (<=0: derived from RANDOM_SEED)");
  CONFIG_ADD_VAR(MATERIAL_INJECT_TAKE, int, 0, "Letters of injected genomes: 0 = taken anywhere (random cells), 1 = own cell/neighbours then anywhere, 2 = own cell/neighbours only");
  CONFIG_ADD_VAR(MATERIAL_AUDIT_INTERVAL, int, 100, "Conservation audit every N updates (0 = off)");
  CONFIG_ADD_VAR(MATERIAL_AUDIT_ABORT, int, 1, "Abort the run on an audit violation");
  CONFIG_ADD_VAR(MATERIAL_CENSUS_FILE, cString, "-", "Prefix of census (.census.gz) and audit (.audit.txt) outputs; - = none");
  CONFIG_ADD_VAR(MATERIAL_CENSUS_INTERVAL, int, 100, "Census row per lineage label every N updates");
  CONFIG_ADD_VAR(MATERIAL_DEBUG_NO_DEATH_RETURN, int, 0, "Negative control: do not return letters on death");
  CONFIG_ADD_VAR(MATERIAL_CENSUS_FOUNDER, cString, "-", "If a sequence is given, census counters of organisms whose genome differs from it go to bucket label+4");
  CONFIG_ADD_VAR(MATERIAL_LIFESPAN, int, 0, "Lifespan constant a in updates since birth (not reset at divide): each organism lives a + U[0,a) updates; 0 = off (use DEATH_METHOD)");
  CONFIG_ADD_VAR(MATERIAL_DEATH_RATE, double, 0.0, "Age-independent death probability per organism per update (material stream)");
  CONFIG_ADD_VAR(COPY_DRAW_RULE, int, 0, "Error draw of h-copy-mat: 0 = draw-first (draw, then seek the letter)\\n1 = find-first (seek the intended letter, draw only if within reach)\\n2 = draw-once (draw first, carry the outcome across stalls)");
  CONFIG_ADD_VAR(COPY_REDRAW_Q, double, -1.0, "Draw-once only: probability that a stall discards the carried draw (<0 = never)");
  CONFIG_ADD_VAR(COPY_MAT_SUB_PROB, double, 0.0, "h-copy-mat: substitution probability per draw (random instruction, uniform)");
  CONFIG_ADD_VAR(COPY_MAT_INS_PROB, double, 0.0, "h-copy-mat: insertion probability per draw");
  CONFIG_ADD_VAR(COPY_MAT_DEL_PROB, double, 0.0, "h-copy-mat: deletion probability per draw");

  CONFIG_ADD_GROUP(DEPRECATED_GROUP,''')

# ---------------- cmake ----------------
edit('../CMakeLists.txt', '  ${MAIN_DIR}/cPopulation.cc\n', '  ${MAIN_DIR}/cMaterial.cc\n  ${MAIN_DIR}/cPopulation.cc\n')

# ---------------- hardware base ----------------
edit('cpu/cHardwareBase.h', '#include <iostream>\n', '#include <iostream>\n#include <vector>\n')
edit('cpu/cHardwareBase.h', '#include "tBuffer.h"\n', '#include "tBuffer.h"\n\nclass cMaterial;\n')
edit('cpu/cHardwareBase.h',
'''  // --------  Organism  ---------
  cOrganism* GetOrganism() { return m_organism; }''',
'''  // --------  Conserved material (h-copy-mat)  ---------
  virtual void MaterialAddHeld(std::vector<long long>& per_kind) const { (void)per_kind; }
  virtual void MaterialReleaseAll(cMaterial& mat, int cell) { (void)mat; (void)cell; }
  virtual void MaterialDiscard() { ; }
  virtual void MaterialClearCarried() { ; }
  virtual int MaterialBucket(cMaterial& mat) { (void)mat; return 0; }
  virtual bool MaterialLifeCheck(cMaterial& mat, int update) { (void)mat; (void)update; return false; }

  // --------  Organism  ---------
  cOrganism* GetOrganism() { return m_organism; }''')

# ---------------- hardware CPU header ----------------
edit('cpu/cHardwareCPU.h',
'  bool Inst_HeadCopy_ifResource(cAvidaContext& ctx);\n',
'''  bool Inst_HeadCopy_ifResource(cAvidaContext& ctx);
  bool Inst_HeadCopyMat(cAvidaContext& ctx);
''')
edit('cpu/cHardwareCPU.h',
'  std::pair<bool, int> m_last_cell_data; //<! If cell data has been previously collected, and it\'s value.\n',
'''  std::pair<bool, int> m_last_cell_data; //<! If cell data has been previously collected, and it's value.

  // conserved material (h-copy-mat)
  std::vector<unsigned char> m_mat_held;   // per memory site: 1 if the site embodies a letter of material
  int m_mat_dpos, m_mat_dk, m_mat_dtype;   // carried draw (draw-once): read position, letter, type
  int m_mat_last_att;                      // time_used at the previous attempt in this gestation (-1 none)
  int m_mat_born, m_mat_life;              // update first seen in the population, drawn lifespan (-1 = not yet)
  int m_mat_founder;                       // genome equals MATERIAL_CENSUS_FOUNDER (1/0, -1 unknown)
  void matSyncHeld();
  bool matDivideAccount(int div_point, int child_size);
public:
  void MaterialAddHeld(std::vector<long long>& per_kind) const;
  void MaterialReleaseAll(cMaterial& mat, int cell);
  void MaterialDiscard() { m_mat_held.clear(); }
  void MaterialClearCarried() { m_mat_dpos = -1; m_mat_last_att = -1; }
  bool MaterialLifeCheck(cMaterial& mat, int update);
  int MaterialBucket(cMaterial& mat);
private:
''')

# ---------------- hardware CPU source ----------------
cc = 'cpu/cHardwareCPU.cc'
edit(cc, '#include "cHardwareCPU.h"\n', '#include "cHardwareCPU.h"\n#include "cMaterial.h"\n')
edit(cc,
'    tInstLibEntry<tMethod>("h-copy-res", &cHardwareCPU::Inst_HeadCopy_ifResource, INST_CLASS_LIFECYCLE, nInstFlag::STALL, "Copy from read-head to write-head if specific resource 1 is available; advance both"),\n',
'''    tInstLibEntry<tMethod>("h-copy-res", &cHardwareCPU::Inst_HeadCopy_ifResource, INST_CLASS_LIFECYCLE, nInstFlag::STALL, "Copy from read-head to write-head if specific resource 1 is available; advance both"),
    tInstLibEntry<tMethod>("h-copy-mat", &cHardwareCPU::Inst_HeadCopyMat, INST_CLASS_LIFECYCLE, 0, "Copy from read-head to write-head taking the letter from conserved material (MATERIAL_MODE); stalls if not within reach"),
''')
edit(cc,
'''  m_memory = *in_seq_p;

  Reset(ctx);                            // Setup the rest of the hardware...
  internalReset();
}''',
'''  m_memory = *in_seq_p;
  m_mat_held.assign(m_memory.GetSize(), 1);
  m_mat_dpos = -1; m_mat_dk = 0; m_mat_dtype = 0; m_mat_last_att = -1; m_mat_born = -1; m_mat_life = -1; m_mat_founder = -1;

  Reset(ctx);                            // Setup the rest of the hardware...
  internalReset();
}''')
edit(cc,
'''  m_mal_active = true;

  return true;
}

int cHardwareCPU::calcCopiedSize''',
'''  m_mal_active = true;
  m_mat_dpos = -1; m_mat_last_att = -1;

  return true;
}

int cHardwareCPU::calcCopiedSize''')
edit(cc,
'''  if (m_world->GetConfig().REQUIRE_EXACT_COPY.Get() && (seq != *offspring_seq) ) {
    return false;
  }

  m_organism->OffspringGenome() = offspring;

  // Cut off everything in this memory past the divide point.
  m_memory.Resize(div_point);

  // Handle Divide Mutations...
  Divide_DoMutations(ctx, mut_multiplier);

  // Many tests will require us to run the offspring through a test CPU;
  // this is, for example, to see if mutations need to be reverted or if
  // lineages need to be updated.
  Divide_TestFitnessMeasures1(ctx);

  if (m_world->GetConfig().DIVIDE_METHOD.Get() != DIVIDE_METHOD_OFFSPRING) {''',
'''  if (m_world->GetConfig().REQUIRE_EXACT_COPY.Get() && (seq != *offspring_seq) ) {
    return false;
  }

  // Conserved material: every kept site must embody a letter; discarded copied sites return.
  const bool mat_on = m_world->GetPopulation().GetMaterial().Enabled() && m_organism->GetCellID() >= 0;
  if (mat_on && !matDivideAccount(div_point, child_size)) return false;

  m_organism->OffspringGenome() = offspring;

  // Cut off everything in this memory past the divide point.
  m_memory.Resize(div_point);
  if (mat_on) m_mat_held.assign(div_point, 1);

  // Handle Divide Mutations...
  Divide_DoMutations(ctx, mut_multiplier);

  // Many tests will require us to run the offspring through a test CPU;
  // this is, for example, to see if mutations need to be reverted or if
  // lineages need to be updated.
  Divide_TestFitnessMeasures1(ctx);

  if (m_world->GetConfig().DIVIDE_METHOD.Get() != DIVIDE_METHOD_OFFSPRING) {''')

impl = r'''

// ======================================================================
//  Conserved material (Code Petri Dish port) — h-copy-mat
// ======================================================================

void cHardwareCPU::matSyncHeld()
{
  const int n = m_memory.GetSize();
  if ((int)m_mat_held.size() < n) m_mat_held.resize(n, 0);          // allocation appends un-copied sites
  else if ((int)m_mat_held.size() > n) {
    m_world->GetPopulation().GetMaterial().Anomaly("held array longer than memory");
    m_mat_held.resize(n);
  }
}

void cHardwareCPU::MaterialAddHeld(std::vector<long long>& per_kind) const
{
  const int n = Apto::Min((int)m_mat_held.size(), m_memory.GetSize());
  for (int i = 0; i < n; i++) if (m_mat_held[i]) per_kind[m_memory[i].GetOp()]++;
}

int cHardwareCPU::MaterialBucket(cMaterial& mat)
{
  int l = m_organism->GetLineageLabel();
  if (l < 0) l = 0;
  if (l > 3) l = 3;
  if (!mat.FounderOn()) return l;
  if (m_mat_founder < 0) m_mat_founder = (cMaterial::KindsOf(m_organism->GetGenome()) == mat.Founder()) ? 1 : 0;
  return m_mat_founder ? l : l + 4;
}

bool cHardwareCPU::MaterialLifeCheck(cMaterial& mat, int update)
{
  if (m_mat_born < 0) { m_mat_born = update; m_mat_life = mat.DrawLifespan(); }
  return update - m_mat_born >= m_mat_life;
}

void cHardwareCPU::MaterialReleaseAll(cMaterial& mat, int cell)
{
  const int n = Apto::Min((int)m_mat_held.size(), m_memory.GetSize());
  for (int i = 0; i < n; i++) if (m_mat_held[i]) mat.Put(cell, m_memory[i].GetOp());
  m_mat_held.clear();
}

// At divide: sites [0, div_point + child_size) are kept (parent and offspring) and must each embody a
// letter; un-copied ones are materialised from reach or the divide is refused. Copied sites past the
// offspring are discarded and returned to the own cell. The offspring's letters go into transit.
bool cHardwareCPU::matDivideAccount(int div_point, int child_size)
{
  cMaterial& mat = m_world->GetPopulation().GetMaterial();
  mat.EnsureInit();
  matSyncHeld();
  const int cell = m_organism->GetCellID();
  cMatCounters& C = mat.Ctr(MaterialBucket(mat));
  const int keep = div_point + child_size;
  std::vector<std::pair<int, int> > taken;
  for (int i = 0; i < keep; i++) {
    if (m_mat_held[i]) continue;
    const int k = m_memory[i].GetOp();
    const int at = mat.Take(cell, k);
    if (at < 0) {
      for (size_t j = 0; j < taken.size(); j++) mat.Put(taken[j].first, taken[j].second);
      C.div_fail_mat++;
      return false;
    }
    taken.push_back(std::make_pair(at, k));
  }
  for (int i = 0; i < keep; i++) m_mat_held[i] = 1;
  C.uncopied_mat += (long long)taken.size();
  for (int i = keep; i < (int)m_mat_held.size(); i++) {
    if (m_mat_held[i]) { mat.Put(cell, m_memory[i].GetOp()); m_mat_held[i] = 0; C.extra_ret++; }
  }
  std::vector<int> kinds;
  kinds.reserve(child_size);
  for (int i = div_point; i < keep; i++) kinds.push_back(m_memory[i].GetOp());
  mat.TransitAdd(kinds);
  C.divides++;
  m_mat_dpos = -1; m_mat_last_att = -1;
  return true;
}

// Mirrors copyOne() of the Code Petri Dish simulator (petri/sim.js v0.3.6, lines 201-227):
//   find-first: seek the intended letter; if not within reach, stall without drawing.
//   draw (own stream): r < pdel -> deletion (no material; read head advances, nothing written)
//                      r < pdel+pins -> insertion of a uniform random letter (read head stays)
//                      r < pdel+pins+psub -> substitution by a uniform random letter
//   draw-once: the outcome (letter, type) is carried across stalls at the same read position;
//              with COPY_REDRAW_Q = q >= 0, each stall discards it with probability q.
//   seek the letter to write (own cell, then neighbours); if none, stall: the cycle is spent,
//   nothing is written, no head moves (under draw-first the next attempt draws again).
bool cHardwareCPU::Inst_HeadCopyMat(cAvidaContext& ctx)
{
  cMaterial& mat = m_world->GetPopulation().GetMaterial();
  const int cell = m_organism->GetCellID();
  if (!mat.Enabled() || cell < 0) return Inst_HeadCopy(ctx);   // material off, or a test CPU
  mat.EnsureInit();

  cHeadCPU& read_head = getHead(nHardware::HEAD_READ);
  cHeadCPU& write_head = getHead(nHardware::HEAD_WRITE);
  read_head.Adjust();
  write_head.Adjust();
  matSyncHeld();

  cMatCounters& C = mat.Ctr(MaterialBucket(mat));
  const int k0 = read_head.GetInst().GetOp();
  const int rh = read_head.GetPosition();
  const int now = m_organism->GetPhenotype().GetTimeUsed();
  C.attempts++;
  if (m_mat_last_att >= 0) { C.gap_sum += now - m_mat_last_att; C.gap_n++; }
  m_mat_last_att = now;

  const int rule = mat.Rule();
  if (rule == 1 && !mat.Has(cell, k0)) { C.stalls++; C.nodraw_stalls++; return true; }

  int k = k0, dtype = 0;
  if (rule == 2 && m_mat_dpos >= 0 && m_mat_dpos == rh) {
    k = m_mat_dk; dtype = m_mat_dtype; C.reused++;
  } else {
    C.draws++;
    const double pd = mat.PDel(), pi = mat.PIns(), ps = mat.PSub();
    const double r = mat.CopyRng().Uniform();
    if (r < pd) {
      C.err_draws++; C.del++;
      ReadInst(k0);
      read_head.Advance();
      m_mat_dpos = -1;
      return true;
    }
    if (r < pd + pi) { k = mat.RandomKind(); dtype = 1; C.err_draws++; }
    else if (r < pd + pi + ps) { k = mat.RandomKind(); dtype = 2; C.err_draws++; C.sub_draws++; }
    if (rule == 2) { m_mat_dk = k; m_mat_dtype = dtype; m_mat_dpos = rh; }
  }

  if (mat.Take(cell, k) < 0) {
    C.stalls++;
    if (dtype != 0) C.err_stalls++;
    if (rule == 2) {
      const double q = mat.Q();
      if (q >= 0 && (q >= 1 || (q > 0 && mat.QRng().Uniform() < q))) m_mat_dpos = -1;
    }
    return true;
  }
  m_mat_dpos = -1;

  const int wp = write_head.GetPosition();
  if (m_mat_held[wp]) { mat.Put(cell, m_memory[wp].GetOp()); C.overwrite_ret++; }
  if (dtype != 1) ReadInst(k0);
  write_head.SetInst(Instruction(k));
  m_mat_held[wp] = 1;
  write_head.SetFlagCopied();
  if (k != k0) { write_head.SetFlagMutated(); write_head.SetFlagCopyMut(); }
  C.written++;
  if (dtype == 2) { C.sub_written++; if (k != k0) C.sub_changed++; }
  if (dtype == 1) { C.ins_written++; write_head.Advance(); }
  else { read_head.Advance(); write_head.Advance(); }
  return true;
}
'''
s = open(cc).read()
s = s.rstrip('\n') + '\n' + impl
open(cc, 'w').write(s)

# ---------------- population ----------------
ph = 'main/cPopulation.h'
edit(ph, 'class cLineage;\n', 'class cLineage;\nclass cMaterial;\n')
edit(ph, '  cBirthChamber birth_chamber;         // Global birth chamber.\n',
     '  cBirthChamber birth_chamber;         // Global birth chamber.\n  cMaterial* m_material;               // Conserved material pools (MATERIAL_MODE)\n')
edit(ph, '  int GetWorldX() const { return world_x; }\n',
     '  int GetWorldX() const { return world_x; }\n  cMaterial& GetMaterial() { return *m_material; }\n')

pc = 'main/cPopulation.cc'
edit(pc, '#include "cPopulation.h"\n', '#include "cPopulation.h"\n#include "cMaterial.h"\n')
edit(pc, ', birth_chamber(world)\n', ', birth_chamber(world)\n, m_material(new cMaterial(world))\n')
edit(pc, '''  delete m_scheduler;
}''', '''  delete m_scheduler;
  delete m_material;
}''')
# death
edit(pc, '''  // And clear it!
  in_cell.RemoveOrganism(ctx);
  if (!organism->IsRunning()) delete organism;''',
'''  // Conserved material: the organism's letters return to its cell.
  if (m_material->Enabled()) m_material->OnDeath(organism, cellID);

  // And clear it!
  in_cell.RemoveOrganism(ctx);
  if (!organism->IsRunning()) delete organism;''')
# placement / drop of offspring
edit(pc, '''      bool org_survived = ActivateOrganism(ctx, offspring_array[i], GetCell(target_cells[i]));''',
'''      if (m_material->Enabled()) m_material->OnPlaced(parent_organism, offspring_array[i]);
      bool org_survived = ActivateOrganism(ctx, offspring_array[i], GetCell(target_cells[i]));''')
edit(pc, '''    } else {
      delete offspring_array[i];
    }
  }
  if (m_world->GetConfig().DIVIDE_METHOD.Get() == DIVIDE_METHOD_SPLIT && parent_alive && m_world->GetConfig().RESET_INPUTS_ON_DIVIDE.Get()) TestForMiniTrace(parent_organism);''',
'''    } else {
      if (m_material->Enabled()) m_material->OnDropped(parent_organism, offspring_array[i], parent_id);
      delete offspring_array[i];
    }
  }
  if (m_world->GetConfig().DIVIDE_METHOD.Get() == DIVIDE_METHOD_SPLIT && parent_alive && m_world->GetConfig().RESET_INPUTS_ON_DIVIDE.Get()) TestForMiniTrace(parent_organism);''')
# inject
edit(pc, '''  // if the injected org already has a group we will assign it to, do not assign group id in activate organism
  if (!inject_group) InjectGenome(cell_id, src, genome, ctx, lineage_label, true);''',
'''  // Conserved material: the injected genome's letters are taken from the pools (or the injection is skipped).
  if (m_material->Enabled() && !m_material->OnInject(ctx, cell_id, genome)) return;

  // if the injected org already has a group we will assign it to, do not assign group id in activate organism
  if (!inject_group) InjectGenome(cell_id, src, genome, ctx, lineage_label, true);''')
# post update
edit(pc, '''  for (int i = 0; i < deme_array.GetSize(); i++) deme_array[i].ProcessUpdate(ctx);
}''', '''  for (int i = 0; i < deme_array.GetSize(); i++) deme_array[i].ProcessUpdate(ctx);

  if (m_material->Enabled()) m_material->PostUpdate(ctx);
}''')

# ---------------- actions ----------------
pa = 'actions/PopulationActions.cc'
edit(pa, '#include "stdlib.h"\n', '#include "stdlib.h"\n#include "cMaterial.h"\n')
actions = r'''

// ---- Conserved material (Code Petri Dish port) ----
class cActionMaterialAudit : public cAction
{
public:
  cActionMaterialAudit(cWorld* world, const cString& args, Feedback&) : cAction(world, args) { ; }
  static const cString GetDescription() { return "Arguments: (none) -- run the conservation audit now"; }
  void Process(cAvidaContext&) { m_world->GetPopulation().GetMaterial().Audit(true); }
};

class cActionMaterialInvade : public cAction
{
private:
  cString m_seq; double m_frac; int m_seed; int m_label;
public:
  cActionMaterialInvade(cWorld* world, const cString& args, Feedback&) : cAction(world, args), m_frac(0.1), m_seed(1), m_label(1)
  {
    cString largs(args);
    m_seq = largs.PopWord();
    if (largs.GetSize()) m_frac = largs.PopWord().AsDouble();
    if (largs.GetSize()) m_seed = largs.PopWord().AsInt();
    if (largs.GetSize()) m_label = largs.PopWord().AsInt();
  }
  static const cString GetDescription() { return "Arguments: <sequence> [frac=0.1] [inj_seed=1] [label=1]"; }
  void Process(cAvidaContext& ctx)
  {
    int n = 0, t = 0, d = 0, s = 0;
    m_world->GetPopulation().GetMaterial().Invade(ctx, (const char*)m_seq, m_frac, (uint64_t)m_seed, m_label, n, t, d, s);
    std::cerr << "[material] invade n=" << n << " target=" << t << " done=" << d << " skipped=" << s << std::endl;
  }
};

class cActionMaterialForkArms : public cAction
{
private:
  cString m_file; int m_parallel;
public:
  cActionMaterialForkArms(cWorld* world, const cString& args, Feedback&) : cAction(world, args), m_parallel(1)
  {
    cString largs(args);
    m_file = largs.PopWord();
    if (largs.GetSize()) m_parallel = largs.PopWord().AsInt();
  }
  static const cString GetDescription() { return "Arguments: <arms_file> [parallel=1]"; }
  void Process(cAvidaContext& ctx) { m_world->GetPopulation().GetMaterial().ForkArms(ctx, (const char*)m_file, m_parallel); }
};
'''
s = open(pa).read()
key = 'void RegisterPopulationActions(cActionLibrary* action_lib)'
i = s.find(key)
if i < 0: sys.exit('RegisterPopulationActions not found')
s = s[:i] + actions.lstrip('\n') + '\n' + s[i:]
s = s.replace('  action_lib->Register<cActionInject>("Inject");\n',
              '  action_lib->Register<cActionInject>("Inject");\n  action_lib->Register<cActionMaterialAudit>("MaterialAudit");\n  action_lib->Register<cActionMaterialInvade>("MaterialInvade");\n  action_lib->Register<cActionMaterialForkArms>("MaterialForkArms");\n', 1)
open(pa, 'w').write(s)
print('patch applied')
