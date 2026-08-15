// Copyright (C) 2010, Guy Barrand. All rights reserved.
// See the file tools.license for terms.

#include "detail/process.hh"
#include "test.hh"
#include <tools/file>
#include <tools/histo/h1d>
#include <tools/mathd>
#include <tools/num2s>

#include <fstream>
#include <random>

#include <tools/rcsv_ntuple>
#include <tools/wcsv_ntuple>

bool test_csv_ntuple(std::ostream& a_out, bool a_verbose)
{
  std::string file;
  if (!tools::num2s(tools::process_id(), file))
  {}
  file = "utest_" + file + ".csv";
  if (a_verbose) a_out << "test_csv_ntuple : file : " << file << std::endl;

  //////////////////////////////////////////////////////////
  unsigned int wentries = 9923;
  if (a_verbose) a_out << "entries : " << wentries << std::endl;

  tools::histo::h1d hgw("rgauss", 100, -5, 5);
  tools::histo::h1d hvdw("rgauss", 100, -5, 5);
  std::vector<double> wsome_entries;

  unsigned int wncols = 0;

  unsigned int iprec = 16;
  double prec = 1e-8;

  ///////////////////////////////////////////////////////////////////////
  /// write : ///////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  if (a_verbose) a_out << "test_csv_ntuple : write ..." << std::endl;

  {
    std::ofstream writer(file.c_str());
    TOOLS_TEST_FUNC(!writer.fail())

    writer.precision(iprec);

    //////////////////////////////////////////////////////////
    /// create and write a ntuple : //////////////////////////
    //////////////////////////////////////////////////////////
    wncols = 0;
    tools::ntuple_booking nbk("rg_rbw", "Randoms");

    nbk.add_column<double>("rgauss");
    wncols++;
    nbk.add_column<float>("rbw");
    wncols++;
    nbk.add_column<std::string>("string");
    wncols++;

    std::vector<float> user_vec_f;
    nbk.add_column<float>("vec_float", user_vec_f);
    wncols++;
    std::vector<double> user_vec_d;
    nbk.add_column<double>("vec_double", user_vec_d);
    wncols++;

    tools::wcsv::ntuple ntuple(writer, a_out, nbk);  // default sep is ','
    ntuple.set_title("csv ntuple");

    TOOLS_TEST_FUNC(ntuple.write_commented_header(a_out))

    if (!tools::equal<unsigned int>(a_out, __FILE__, __LINE__, ntuple.columns().size(), wncols))
      return false;

    tools::wcsv::ntuple::column<double>* col_rgauss = ntuple.find_column<double>("rgauss");
    if (!tools::valid_pointer(a_out, __FILE__, __LINE__, col_rgauss)) return false;
    tools::wcsv::ntuple::column<float>* col_rbw = ntuple.find_column<float>("rbw");
    if (!tools::valid_pointer(a_out, __FILE__, __LINE__, col_rbw)) return false;
    tools::wcsv::ntuple::column<std::string>* col_str = ntuple.find_column<std::string>("string");
    if (!tools::valid_pointer(a_out, __FILE__, __LINE__, col_str)) return false;

    std::mt19937 gen;
    std::normal_distribution<double> vrg(5.0, 2.0);
    std::normal_distribution<double> rg(1.0, 2.0);
    std::cauchy_distribution<float> rbwf(0.0, 0.5);

    std::string stmp;
    for (unsigned int count = 0; count < wentries; count++)
    {
      double vd = rg(gen);
      hgw.fill(vd);
      if (!tools::equal(a_out, __FILE__, __LINE__, col_rgauss->fill(vd), true)) return false;
      if (!tools::equal(a_out, __FILE__, __LINE__, col_rbw->fill(rbwf(gen)), true)) return false;
      if (!tools::num2s(count, stmp))
      {}
      if (!tools::equal(a_out, __FILE__, __LINE__, col_str->fill("str " + stmp), true))
        return false;

      {
        double dnumber = vrg(gen);
        size_t number = size_t(dnumber >= 0 ? dnumber : 0);
        user_vec_f.resize(number);
        for (size_t i = 0; i < number; i++)
          user_vec_f[i] = rbwf(gen);
      }

      {
        double dnumber = vrg(gen);
        size_t number = size_t(dnumber >= 0 ? dnumber : 0);
        user_vec_d.resize(number);
        for (size_t i = 0; i < number; i++)
        {
          double v = rg(gen);
          user_vec_d[i] = v;
          hvdw.fill(v);
        }
      }

      if (!tools::equal(a_out, __FILE__, __LINE__, ntuple.add_row(), true)) return false;

      if (count <= 5)
      {
        if (a_verbose) a_out << "count " << count << " : " << vd << std::endl;
        wsome_entries.push_back(vd);
      }
    }
  }

  ///////////////////////////////////////////////////////////////////////
  /// read : ////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////

  if (a_verbose) a_out << "test_csv_ntuple : read ..." << std::endl;

  {
    std::ifstream reader(file.c_str());
    TOOLS_TEST_FUNC(!reader.fail())

    reader.precision(iprec);

    tools::rcsv::ntuple ntuple(reader);

    TOOLS_TEST_FUNC(ntuple.initialize_from_commented_header(a_out))

    if (!tools::equal<unsigned int>(a_out, __FILE__, __LINE__, ntuple.columns().size(), wncols))
      return false;

    tools::read::icolumn<double>* col_rgauss = ntuple.find_column<double>("rgauss");
    if (!tools::valid_pointer(a_out, __FILE__, __LINE__, col_rgauss)) return false;

    tools::read::icolumn<float>* col_rbw = ntuple.find_column<float>("rbw");
    if (!tools::valid_pointer(a_out, __FILE__, __LINE__, col_rbw)) return false;

    tools::read::icolumn<std::string>* col_str = ntuple.find_column<std::string>("string");
    if (!tools::valid_pointer(a_out, __FILE__, __LINE__, col_str)) return false;

    tools::read::icolumn<std::vector<float>>* col_vec_float =
      ntuple.find_column<std::vector<float>>("vec_float");
    if (!tools::valid_pointer(a_out, __FILE__, __LINE__, col_vec_float)) return false;

    tools::read::icolumn<std::vector<double>>* col_vec_double =
      ntuple.find_column<std::vector<double>>("vec_double");
    if (!tools::valid_pointer(a_out, __FILE__, __LINE__, col_vec_double)) return false;

    std::vector<double> user_vec_d;
    std::string stmp;
    tools::histo::h1d hgr("rgauss", 100, -5, 5);
    tools::histo::h1d hvdr("rgauss", 100, -5, 5);
    std::vector<double> rsome_entries;

    ntuple.start();
    tools::uint64 row = 0;
    while (ntuple.next())
    {
      double vd;
      TOOLS_TEST_FUNC(col_rgauss->get_entry(vd))
      hgr.fill(vd);

      float vf;
      TOOLS_TEST_FUNC(col_rbw->get_entry(vf))

      std::string vs;
      TOOLS_TEST_FUNC(col_str->get_entry(vs))
      stmp = "str ";
      if (!tools::numas(row, stmp))
      {}
      if (!tools::equal(a_out, __FILE__, __LINE__, vs, stmp)) return false;

      if (!tools::equal(a_out, __FILE__, __LINE__, col_vec_double->get_entry(user_vec_d), true))
        return false;
      {
        tools_vforcit(double, user_vec_d, it) hvdr.fill(*it);
      }

      if (row <= 5)
      {
        if (a_verbose) a_out << "count " << row << " : " << vd << std::endl;
        rsome_entries.push_back(vd);
      }
      row++;
    }

    tools::uint64 entries = row;
    if (!tools::equal<tools::uint64>(a_out, __FILE__, __LINE__, entries, tools::uint64(wentries)))
      return false;

    if (!tools::equal(a_out, __FILE__, __LINE__, hgw.equals(hgr, prec, ::fabs), true))
    {
      hgw.hprint(a_out);
      hgr.hprint(a_out);
      return false;
    }
    if (!tools::equal(a_out, __FILE__, __LINE__,
                      tools::vectors_are_equal(wsome_entries, rsome_entries, prec, tools::dfabs),
                      true))
      return false;

    if (!tools::equal(a_out, __FILE__, __LINE__, hvdw.equals(hvdr, prec, ::fabs), true))
      return false;
  }

  ///////////////////////////////////////////////////////////////////////
  TOOLS_TEST_FUNC(tools::file::std_remove(file))
  return true;
}
