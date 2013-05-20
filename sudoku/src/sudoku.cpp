/**
 * @file sudoku.cpp
 *
 * @brief sophisticated sudoku solver logic.
 *
 * @author Andreas Lobbes
 *         | alobbes@web.de
 *         | andreas.lobbes@gmail.com
 *         | http://andreas-lobbes.blogspot.com
 *
 * @date $Date: 2013/04/02 11:39:35 $
 * @version $Revision: 1.33 $
 *
 * @copyright (C) 2013 Andreas Lobbes
 *
 * @see
 * - http://en.wikipedia.org/wiki/Sudoku
 * - http://en.wikipedia.org/wiki/Mathematics_of_Sudoku
 * - http://school.maths.uwa.edu.au/~gordon/sudokumin
 */
#include <iostream>
#include <sstream>
#include <iosfwd>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <cctype>
#include <set>
#include <sys/types.h>

#ifdef _MSC_VER
#define __int64_t __int64
#endif


namespace Sudoku {

class Cell;
class Square;
class Grid;

class Cell {
  friend class Grid;

  /**
   * @brief cell constructor.
   *
   * @param id      cell id.
   * @param rowmap  map of fixed numbers this row.
   * @param colmap  map of fixed numbers this column.
   * @param sqrmap  map of fixed numbers this square.
   * @param srwmap  map of located (but not fixed) numbers this row of this square.
   * @param sclmap  map of located (but not fixed) numbers this column of this square.
   * @param qr0map  map of located (but not fixed) numbers this row of neighboring square.
   * @param qr1map  map of located (but not fixed) numbers this row of (2nd) neighboring square.
   * @param qc0map  map of located (but not fixed) numbers this column of neighboring square.
   * @param qc1map  map of located (but not fixed) numbers this column of (2nd) neighboring square.
   * @param sr0map  map of located (but not fixed) numbers of neighboring row of this square.
   * @param sr1map  map of located (but not fixed) numbers of (2nd) neighboring row of this square.
   * @param sc0map  map of located (but not fixed) numbers of neighboring column of this square.
   * @param sc1map  map of located (but not fixed) numbers of (2nd) neighboring column of this square.
   * @param _xr0map cache map for auto-fixing numbers this row. @see Sudoku::Grid::_xrmap
   * @param _xr1map cache map for auto-fixing numbers this row. @see Sudoku::Grid::_xrmap
   * @param _xc0map cache map for auto-fixing numbers this column. @see Sudoku::Grid::_xcmap
   * @param _xc1map cache map for auto-fixing numbers this column. @see Sudoku::Grid::_xcmap
   */
  Cell(unsigned id, unsigned &rowmap, unsigned &colmap, unsigned &sqrmap,
      unsigned &srwmap, unsigned &sclmap,
      unsigned &qr0map, unsigned &qr1map, unsigned &qc0map, unsigned &qc1map,
      unsigned &sr0map, unsigned &sr1map, unsigned &sc0map, unsigned &sc1map,
      unsigned &_xr0map, unsigned &_xr1map, unsigned &_xc0map, unsigned &_xc1map):
      next(nullptr), id(id),
      rowmap(rowmap), colmap(colmap), sqrmap(sqrmap), fixbit(0),
      srwmap(srwmap), sclmap(sclmap),
      qr0map(qr0map), qr1map(qr1map), qc0map(qc0map), qc1map(qc1map),
      sr0map(sr0map), sr1map(sr1map), sc0map(sc0map), sc1map(sc1map),
      _xr0map(_xr0map), _xr1map(_xr1map), _xc0map(_xc0map), _xc1map(_xc1map) {
  }

public:
  bool IsValid() const {
    return IsFixed() || ((rowmap | colmap | sqrmap | qr0map | qr1map | qc0map | qc1map) != ((1 << 9) - 1));
  }
  bool IsFixed() const {
    return fixbit != 0;
  }
  unsigned GetFixbit() const {
    return fixbit;
  }
  bool CanFix(unsigned bit) const {
    return !((rowmap | colmap | sqrmap | qr0map | qr1map | qc0map | qc1map) & bit);
  }
  int AutoFix() {
    if (IsFixed())
      return 0;
    switch (rowmap | colmap | sqrmap | qr0map | qr1map | qc0map | qc1map) {
    case 0x1fe: return (int) Fix(0x01);
    case 0x1fd: return (int) Fix(0x02);
    case 0x1fb: return (int) Fix(0x04);
    case 0x1f7: return (int) Fix(0x08);
    case 0x1ef: return (int) Fix(0x10);
    case 0x1df: return (int) Fix(0x20);
    case 0x1bf: return (int) Fix(0x40);
    case 0x17f: return (int) Fix(0x80);
    case 0xff:  return (int) Fix(0x100);
    case 0x1ff: return -1;
    }
    unsigned bit;
    switch(bit = (srwmap & sclmap)) {
    case 0:
      return 0;
    case 0x01:
    case 0x02:
    case 0x04:
    case 0x08:
    case 0x10:
    case 0x20:
    case 0x40:
    case 0x80:
    case 0x100:
      return CanFix(bit) ? (int) Fix(bit) : -1;
    default:
      return -1;
    }
  }
  unsigned Fix(unsigned bit) {
    rowmap |= bit;
    colmap |= bit;
    sqrmap |= bit;
    return fixbit = bit;
  }
  unsigned Unfix() {
    unsigned bit = fixbit;
    unsigned mask = ~bit;
    rowmap &= mask;
    colmap &= mask;
    sqrmap &= mask;
    fixbit = 0;
    return bit;
  }

private:
  Cell *next;           /**< @brief used for chaining. */
  unsigned id;          /**< @brief cell id. */
  unsigned &rowmap;     /**< @brief map of fixed numbers this row. */
  unsigned &colmap;     /**< @brief map of fixed numbers this column. */
  unsigned &sqrmap;     /**< @brief map of fixed numbers this square. */
  unsigned fixbit;      /**< @brief fixed number of this cell (1 << (n - 1)). */
  unsigned &srwmap;     /**< @brief map of located (but not fixed) numbers this row of this square. */
  unsigned &sclmap;     /**< @brief map of located (but not fixed) numbers this column of this square. */
  unsigned &qr0map;     /**< @brief map of located (but not fixed) numbers this row of neighboring square. */
  unsigned &qr1map;     /**< @brief map of located (but not fixed) numbers this row of (2nd) neighboring square. */
  unsigned &qc0map;     /**< @brief map of located (but not fixed) numbers this column of neighboring square. */
  unsigned &qc1map;     /**< @brief map of located (but not fixed) numbers this column of (2nd) neighboring square. */
  unsigned &sr0map;     /**< @brief map of located (but not fixed) numbers of neighboring row of this square. */
  unsigned &sr1map;     /**< @brief map of located (but not fixed) numbers of (2nd) neighboring row of this square. */
  unsigned &sc0map;     /**< @brief map of located (but not fixed) numbers of neighboring column of this square. */
  unsigned &sc1map;     /**< @brief map of located (but not fixed) numbers of (2nd) neighboring column of this square. */
  unsigned &_xr0map;    /**< @brief cache map for auto-fixing numbers this row. @see Sudoku::Grid::_xrmap */
  unsigned &_xr1map;    /**< @brief cache map for auto-fixing numbers this row. @see Sudoku::Grid::_xrmap */
  unsigned &_xc0map;    /**< @brief cache map for auto-fixing numbers this column. @see Sudoku::Grid::_xcmap */
  unsigned &_xc1map;    /**< @brief cache map for auto-fixing numbers this column. @see Sudoku::Grid::_xcmap */
};

class Square {
  friend class Grid;

  Square(unsigned &r0map, unsigned &r1map, unsigned &r2map,
        unsigned &c0map, unsigned &c1map, unsigned &c2map,
        unsigned &sqrmap,
        unsigned &sr0map, unsigned &sr1map, unsigned &sr2map,
        unsigned &sc0map, unsigned &sc1map, unsigned &sc2map
        ): sqrmap(sqrmap) {
    rowmap[0] = &r0map;
    rowmap[1] = &r1map;
    rowmap[2] = &r2map;
    colmap[0] = &c0map;
    colmap[1] = &c1map;
    colmap[2] = &c2map;

    srwmap[0] = &sr0map;
    srwmap[1] = &sr1map;
    srwmap[2] = &sr2map;
    sclmap[0] = &sc0map;
    sclmap[1] = &sc1map;
    sclmap[2] = &sc2map;
  };

  unsigned *rowmap[3];  /**< @brief map of fixed numbers of rows crossing this square */
  unsigned *colmap[3];  /**< @brief map of fixed numbers of columns crossing this square */
  unsigned &sqrmap;     /**< @brief map of fixed numbers this square */
  unsigned *srwmap[3];  /**< @brief map of located (but not fixed) numbers of rows this square */
  unsigned *sclmap[3];  /**< @brief map of located (bit not fixed) numbers of columns this square */
  unsigned _srmap[3];   /**< @brief cache map of new located (but not fixed) numbers of rows this square */
  unsigned _scmap[3];   /**< @brief cache map of new located (bit not fixed) numbers of columns this square */
};

class Grid {
public:
  class Data {
    friend std::ostream &operator<<(std::ostream &os, const Data &gd);
    friend std::istream &operator>>(std::istream &is, Data &gd);
  public:
    Data() {
      std::memset(cells, 0, sizeof(cells));
    }
    Data(const Data &gd) {
      *this = gd;
    }
    Data(const Grid &g) {
      *this = g;
    }
    Data &operator=(const Data &gd) {
      std::memmove(cells, gd.cells, sizeof(cells));
      return *this;
    }
    Data &operator=(const Grid &g) {
      for (int i = 0; i < 81; ++i) {
        Cell *cell = g[i];
        char c;
        switch(cell->GetFixbit()) {
        default:
        case 0:         c = 0; break;
        case 1 << 0:    c = 1; break;
        case 1 << 1:    c = 2; break;
        case 1 << 2:    c = 3; break;
        case 1 << 3:    c = 4; break;
        case 1 << 4:    c = 5; break;
        case 1 << 5:    c = 6; break;
        case 1 << 6:    c = 7; break;
        case 1 << 7:    c = 8; break;
        case 1 << 8:    c = 9; break;
        }
        cells[i] = c;
      }
      return *this;
    }

    bool operator==(const Data &gd) const {
      return std::memcmp(cells, gd.cells, sizeof(cells)) == 0;
    }
    bool operator!=(const Data &gd) const {
      return !(*this == gd);
    }
    bool operator<(const Data &gd) const {
      return std::memcmp(cells, gd.cells, sizeof(cells)) < 0;
    }
    bool operator>(const Data &gd) const {
      return std::memcmp(cells, gd.cells, sizeof(cells)) > 0;
    }

    int operator[](int idx) const {
      return cells[idx];
    }

    Data &Transpose(Data &gd, int n) const {
      if (!n)
        return gd = *this;

      for (int i = 0; i < 9; ++i) {
        for (int j = 0; j < 9; ++j) {
          gd.grid[j][i] = grid[i][j];
        }
      }
      return gd;
    }
    Data &QRowOrder(Data &gd, int n) const {
      unsigned perm = Perm3(n);
      for (int qr = 0; qr < 9; perm >>= 8, qr += 3) {
        const char *q = &(grid[qr][0]);
        char *p = &(gd.grid[((char) perm) * 3][0]);
        for (int i = 27; --i >= 0;) {
          *p++ = *q++;
        }
      }
      return gd;
    }
    Data &QColOrder(Data &gd, int n) const {
      unsigned perm = Perm3(n);
      for (int qc = 0; qc < 9; perm >>= 8, qc += 3) {
        const char *q = &(grid[0][qc]);
        char *p = &(gd.grid[0][((char) perm) * 3]);
        for (int r = 9; --r >= 0; p += 6, q += 6) {
          *p++ = *q++;
          *p++ = *q++;
          *p++ = *q++;
        }
      }
      return gd;
    }
    Data &RowOrder(Data &gd, int n0, int n1, int n2) const {
      unsigned perm = Perm3(n0);
      const char *q = &(grid[0][0]);
      for (int r = 3; --r >= 0; perm >>= 8) {
        char *p = &(gd.grid[(int)((char) perm)][0]);
        for (int c = 9; --c >= 0;)
          *p++ = *q++;
      }
      perm = Perm3(n1);
      q = &(grid[3][0]);
      for (int r = 3; --r >= 0; perm >>= 8) {
        char *p = &(gd.grid[3 + (char) perm][0]);
        for (int c = 9; --c >= 0;)
          *p++ = *q++;
      }
      perm = Perm3(n2);
      q = &(grid[6][0]);
      for (int r = 3; --r >= 0; perm >>= 8) {
        char *p = &(gd.grid[6 + (char) perm][0]);
        for (int c = 9; --c >= 0;)
          *p++ = *q++;
      }
      return gd;
    }
    Data &ColOrder(Data &gd, int n0, int n1, int n2) const {
      unsigned perm = Perm3(n0);
      for (int c = 0; c < 3; perm >>= 8, ++c) {
        const char *q = &(grid[0][c]);
        char *p = &(gd.grid[0][(int)((char) perm)]);
        for (int r = 9; --r >= 0; p += 9, q += 9) {
          *p = *q;
        }
      }
      perm = Perm3(n1);
      for (int c = 0; c < 3; perm >>= 8, ++c) {
        const char *q = &(grid[0][3 + c]);
        char *p = &(gd.grid[0][3 + (char) perm]);
        for (int r = 9; --r >= 0; p += 9, q += 9) {
          *p = *q;
        }
      }
      perm = Perm3(n2);
      for (int c = 0; c < 3; perm >>= 8, ++c) {
        const char *q = &(grid[0][6 + c]);
        char *p = &(gd.grid[0][6 + (char) perm]);
        for (int r = 9; --r >= 0; p += 9, q += 9) {
          *p = *q;
        }
      }
      return gd;
    }
    Data &LowOrder(Data &gd) const {
      char perm[9];
      unsigned map = 0;
      char mapc = 0;
      const char *q = cells;
      char *p = gd.cells;
      for (int i = 81; --i >= 0;) {
        char c = *q++;
        if (c != 0) {
          unsigned m = 1U << (c - 1);
          if (!(map & m)) {
            map |= m;
            perm[(int) c] = ++mapc;
          }
          c = perm[(int) c];
        }
        *p++ = c;
      }
      return gd;
    }
    Data &Normalize(Data &gd) const {
      Data T, QR, QC, R, C, L;
      gd = *this;
      for (int i = 0; i < 2; ++i) {
        Transpose(T, i);
        for (int j = 0; j < 6; ++j) {
          T.QRowOrder(QR, j);
          for (int k = 0; k < 6; ++k) {
            QR.QColOrder(QC, k);
            for (int r0 = 0; r0 < 6; ++r0) {
              for (int r1 = 0; r1 < 6; ++r1) {
                for (int r2 = 0; r2 < 6; ++r2) {
                  QC.RowOrder(R, r0, r1, r2);
                  for (int c0 = 0; c0 < 6; ++c0) {
                    for (int c1 = 0; c1 < 6; ++c1) {
                      for (int c2 = 0; c2 < 6; ++c2) {
                        R.ColOrder(C, c0, c1, c2);
                        C.LowOrder(L);
                        if (L < gd)
                          gd = L;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      return gd;
    }

    std::ostream &Print(std::ostream &os) const {
      for (int row = 0; row < 9; ++row) {
        if ((row != 0) && ((row % 3) == 0))
          os << std::endl;
        for (int col = 0; col < 9; ++col) {
          if (col != 0) {
            os << " ";
            if ((col % 3) == 0)
              os << " ";
          }
          char c;
          if ((c = grid[row][col]) != 0)
            c += '0';
          else
            c = '.';
          os << c;
        }
        os << std::endl;
      }
      return os;
    }

  private:
    static unsigned Perm3(int n) {      // g++ still does not support const init.
      switch(n) {
      default:
      case 0:   return 0x020100U;
      case 1:   return 0x010200U;
      case 2:   return 0x020001U;
      case 3:   return 0x000201U;
      case 4:   return 0x010002U;
      case 5:   return 0x000102U;
      }
    }

    union {
      char grid[9][9];
      char cells[81];
    };
  };

  Grid() {
    Init();
  }
  Grid(const Grid &src) {
    Init();
    *this = src;
  }
  Grid(const Grid::Data &src) {
    Init();
    *this = src;
  }
  ~Grid() {
    for (int row = 0; row < 9; ++row) {
      for (int col = 0; col < 9; ++col) {
        delete grid[row][col];
      }
    }
    for (int qrow = 0; qrow < 3; ++qrow) {
      for (int qcol = 0; qcol < 3; ++qcol) {
        delete sgrid[qrow][qcol];
      }
    }
  }

  bool IsValid() {
    Grid g;
    for (int i = 0; i < 81; ++i) {
      unsigned b = cells[i]->GetFixbit();
      Cell *c = g.cells[i];
      if (!c->CanFix(b))
        return false;
      c->Fix(b);
    }
    __int64_t N;
    Permutate(N, [](__int64_t, const Grid&) {return false;});
    return N != 0;
  }
  bool IsUnique() {
    __int64_t N;
    Permutate(N, [](__int64_t N, const Grid&) {return N == 1;});
    return N == 1;
  }
  explicit operator bool() {
    return IsValid();
  }

  Grid &Reset() {
    std::memset(rowmap, 0, sizeof(rowmap));
    std::memset(colmap, 0, sizeof(colmap));
    std::memset(sqrmap, 0, sizeof(sqrmap));
    std::memset(qrwmap, 0, sizeof(qrwmap));
    std::memset(qclmap, 0, sizeof(qclmap));
    for (int i = 0; i < 81; ++i)
      cells[i]->fixbit = 0;
    return *this;
  }
  Grid &Random() {
    std::srand((unsigned) std::time(nullptr));
    Reset();
    unsigned map = 0;
    for (int qr = 0; qr < 3; ++qr) {
      for (int qc = 0; qc < 3; ++qc) {
        int r = std::rand() % 3;
        int c = std::rand() % 3;
        unsigned b;
        do b = 1U << (std::rand() % 9);
        while (map & b);
        map |= b;
        grid[qr * 3 + r][qc * 3 + c]->Fix(b);
      }
    }
    return *this;
  }
  Grid &operator=(const Grid &src) {
    memcpy(rowmap, src.rowmap, sizeof(rowmap));
    memcpy(colmap, src.colmap, sizeof(colmap));
    memcpy(sqrmap, src.sqrmap, sizeof(sqrmap));
    memcpy(qrwmap, src.qrwmap, sizeof(qrwmap));
    memcpy(qclmap, src.qclmap, sizeof(qclmap));
    for (int i = 0; i < 81; ++i)
      cells[i]->fixbit = src.cells[i]->fixbit;
    return *this;
  }
  Grid &operator=(const Data &src) {
    Reset();
    for (int i = 0; i < 81; ++i) {
      int b = src[i];
      if (b)
        cells[i]->Fix(1U << (b - 1));
    }
    return *this;
  }
  Cell *operator[](int idx) const {
    return cells[idx];
  }

  template<typename EnumF>
  bool Permutate(__int64_t &N, EnumF ef, bool r0 = true) {
    if (r0) {
      ClearStk();
      N = 0;
    }

    Stkp stkp;
    if (AutoFix(stkp) < 0)
      return true;

    Cell *cell;
    if ((cell = cunfix) != nullptr) {
      cunfix = cell->next;
      cell->next = cfixed;
      cfixed = cell;
      for (unsigned bit = 1; !(bit & (1 << 9)); bit <<= 1) {
        if (!cell->CanFix(bit))
          continue;
        cell->Fix(bit);
        bool go = Permutate(N, ef, false);
        cell->Unfix();
        if (!go) {
          Unfix(stkp);
          return false;
        }
      }
      Unfix(stkp);
      return true;
    }

    ++N;
    bool go = ef(N, *this);
    Unfix(stkp);
    return go;
  }
  bool Permutate(__int64_t &N) {
    return Permutate(N, [](__int64_t, const Grid&) {return true;});
  }

  std::ostream &Print(std::ostream &os) const {
    return Data(*this).Print(os);
  }

private:
  void Init() {
    memset(cells, 0, sizeof(cells));
    for (int row = 0; row < 9; ++row) {
      for (int col = 0; col < 9; ++col) {
        grid[row][col] = new Cell(row * 9 + col,
            rowmap[row], colmap[col], sqrmap[row / 3][col / 3],
            qrwmap[row][col / 3],
            qclmap[col][row / 3],
            qrwmap[row][((col / 3) + 1) % 3],
            qrwmap[row][((col / 3) + 2) % 3],
            qclmap[col][((row / 3) + 1) % 3],
            qclmap[col][((row / 3) + 2) % 3],
            qrwmap[(row / 3) * 3 + (row + 1) % 3][col / 3],
            qrwmap[(row / 3) * 3 + (row + 2) % 3][col / 3],
            qclmap[(col / 3) * 3 + (col + 1) % 3][row / 3],
            qclmap[(col / 3) * 3 + (col + 2) % 3][row / 3],
            _xrmap[row][0], _xrmap[row][1], _xcmap[col][0], _xcmap[col][1]
        );
      }
    }
    for (int qrow = 0; qrow < 3; ++qrow) {
      for (int qcol = 0; qcol < 3; ++qcol) {
        sgrid[qrow][qcol] = new Square(
            rowmap[3 * qrow], rowmap[3 * qrow + 1], rowmap[3 * qrow + 2],
            colmap[3 * qcol], colmap[3 * qcol + 1], colmap[3 * qcol + 2],
            sqrmap[qrow][qcol],
            qrwmap[3 * qrow][qcol], qrwmap[3 * qrow + 1][qcol], qrwmap[3 * qrow + 2][qcol],
            qclmap[3 * qcol][qrow], qclmap[3 * qcol + 1][qrow], qclmap[3 * qcol + 2][qrow]);
      }
    }
    Reset();
  }

  struct Stkp {
    Cell *cfixed;
    unsigned qrwmap[9][3];
    unsigned qclmap[9][3];
  };

  void ClearStk() {
    cfixed = nullptr;
    cunfix = nullptr;
    std::memset(qrwmap, 0, sizeof(qrwmap));
    std::memset(qclmap, 0, sizeof(qclmap));
    for (unsigned i = 0; i < 81; ++i) {
      Cell *cell = cells[i];
      Cell *&chain = cell->IsFixed() ? cfixed : cunfix;
      cell->next = chain;
      chain = cell;
    }
  }
  int AutoFix(Stkp &stkp) {
    stkp.cfixed = cfixed;
    memcpy(stkp.qrwmap, qrwmap, sizeof(qrwmap));
    memcpy(stkp.qclmap, qclmap, sizeof(qclmap));
    int m, n = 0;
    do {
      // locate unfixed numbers within rows / columns of a square
      for (int i = 0; i < 9; ++i) {
        qrwmap[i][0] = qrwmap[i][1] = qrwmap[i][2] = 0x1ff;
        qclmap[i][0] = qclmap[i][1] = qclmap[i][2] = 0x1ff;
      }
      Cell *cell, **it = &cfixed;
      while ((cell = *it) != nullptr) {
        it = &cell->next;
        unsigned xmap = ~cell->fixbit;
        cell->sr0map &= xmap;
        cell->sr1map &= xmap;
        cell->sc0map &= xmap;
        cell->sc1map &= xmap;
      }
      it = &cunfix;
      while ((cell = *it) != nullptr) {
        it = &cell->next;
        unsigned xmap = cell->rowmap | cell->colmap;
        cell->sr0map &= xmap;
        cell->sr1map &= xmap;
        cell->sc0map &= xmap;
        cell->sc1map &= xmap;
      }
      // locate unfixed numbers within rows / columns.
      // the less cell numbers are fixed already, the more expensive is this logic.
      // but nevertheless this logic boosts overall performance when solving sudokus.
      for (int i = 0; i < 9; ++i) {
        _xrmap[i][0] = _xcmap[i][0] = 0;
        _xrmap[i][1] = _xcmap[i][1] = 0x1ff;
      }
      it = &cunfix;
      while ((cell = *it) != nullptr) {
        it = &cell->next;
        unsigned xmap = ~(cell->rowmap | cell->colmap | cell->sqrmap);
        unsigned ymap = (cell->_xr0map ^ xmap) & xmap;
        cell->_xr0map |= ymap;
        cell->_xr1map &= ~(xmap ^ ymap);
        ymap = (cell->_xc0map ^ xmap) & xmap;
        cell->_xc0map |= ymap;
        cell->_xc1map &= ~(xmap ^ ymap);
      }
      it = &cunfix;
      while ((cell = *it) != nullptr) {
        it = &cell->next;
        unsigned xmap = ~(cell->rowmap | cell->colmap | cell->sqrmap);
        unsigned bit = (cell->_xr1map | cell->_xc1map) & xmap;
        if (bit == 0)
          continue;
        switch(bit) {
        default:
          Unfix(stkp);
          return -1;
        case 0x01:
        case 0x02:
        case 0x04:
        case 0x08:
        case 0x10:
        case 0x20:
        case 0x40:
        case 0x80:
        case 0x100:
          cell->srwmap |= bit;
          cell->sclmap |= bit;
          break;
        }
      }
      // more expensive logic to locate further unfixed numbers within rows / columns of a square.
      // effectively this slightly decreases overall performance when solving sudokus.
      for (int qr = 0; qr < 3; ++qr) {
        for (int qc = 0; qc < 3; ++qc) {
          Square *s = sgrid[qr][qc];
          for (int i = 0; i < 3; ++i) {
            Square *sr0 = sgrid[qr][(qc + 1) % 3];
            Square *sr1 = sgrid[qr][(qc + 2) % 3];
            Square *sc0 = sgrid[(qr + 1) % 3][qc];
            Square *sc1 = sgrid[(qr + 2) % 3][qc];
            int i0 = (i + 1) % 3;
            int i1 = (i + 2) % 3;
            s->_srmap[i] =
                ((*s->rowmap[i0] | *sr0->srwmap[i0]) &  (*s->rowmap[i1] | *sr1->srwmap[i1]) & (~s->sqrmap)) |
                ((*s->rowmap[i1] | *sr0->srwmap[i1]) &  (*s->rowmap[i0] | *sr1->srwmap[i0]) & (~s->sqrmap));
            s->_scmap[i] =
                ((*s->colmap[i0] | *sc0->sclmap[i0]) &  (*s->colmap[i1] | *sc1->sclmap[i1]) & (~s->sqrmap)) |
                ((*s->colmap[i1] | *sc0->sclmap[i1]) &  (*s->colmap[i0] | *sc1->sclmap[i0]) & (~s->sqrmap));
          }
        }
      }
      for (int qr = 0; qr < 3; ++qr) {
        for (int qc = 0; qc < 3; ++qc) {
          Square *s = sgrid[qr][qc];
          for (int i = 0; i < 3; ++i) {
            *s->srwmap[i] |= s->_srmap[i];
            *s->sclmap[i] |= s->_scmap[i];
          }
        }
      }

      // automatically fix cell's numbers where they can be derived.
      m = 0;
      it = &cunfix;
      while ((cell = *it) != nullptr) {
        int bit = cell->AutoFix();
        if (bit < 0) {
          Unfix(stkp);
          return -1;
        }
        if (bit == 0) {
          it = &cell->next;
          continue;
        }
        *it = cell->next;
        cell->next = cfixed;
        cfixed = cell;
        ++m;
      }
      n += m;
    } while (m != 0);
    return n;
  }
  void Unfix(Stkp &stkp) {
    Cell *cell;
    while ((cell = cfixed) != stkp.cfixed) {
      cfixed = cell->next;
      cell->Unfix();
      cell->next = cunfix;
      cunfix = cell;
    }
    memcpy(qrwmap, stkp.qrwmap, sizeof(qrwmap));
    memcpy(qclmap, stkp.qclmap, sizeof(qclmap));
  }

  union {
    Cell *cells[81];
    Cell *grid[9][9];
  };
  Square *sgrid[3][3];
  Cell *cfixed;                 /**< @brief chain of already fixed cells. */
  Cell *cunfix;                 /**< @brief chain of still unfixed cells. */
  unsigned rowmap[9];           /**< @brief map of fixed numbers per row. */
  unsigned colmap[9];           /**< @brief map of fixed numbers per column. */
  unsigned sqrmap[3][3];        /**< @brief map of fixed numbers per square. */
  unsigned qrwmap[9][3];        /**< @brief map of located (but not fixed) numbers per row and square. */
  unsigned qclmap[9][3];        /**< @brief map of located (but not fixed) numbers per column and square */

  /**
   * cache map to auto-fix row numbers.
   *
   * _xrmap[row][0] accumulates all allocatable bits.
   * at end of iteration (of all row cells) this becomes effectively
   * _xrmap[row][0] = ~rowmap[row] & 0x1ff.
   *
   * _xrmap[row][1] filters all bits, they were allocatable once only
   * by one cell of iterated row.
   */
  unsigned _xrmap[9][2];
  /**
   * cache map to auto-fix column numbers.
   *
   * _xcmap[col][0] accumulates all allocatable bits.
   * at end of iteration (of all column cells) this becomes effectively
   * _xcmap[col][0] = ~colmap[col] & 0x1ff.
   *
   * _xcmap[col][1] filters all bits, they were allocatable once only
   * by one cell of iterated column.
   */
  unsigned _xcmap[9][2];
};

std::ostream &operator<<(std::ostream &os, const Sudoku::Grid &grid) {
  return os << Sudoku::Grid::Data(grid);
}

std::istream &operator>>(std::istream &is, Sudoku::Grid &grid) {
  Sudoku::Grid::Data gd;
  if (!(is >> gd) || !Sudoku::Grid(gd)) {
    is.setstate(std::ios::failbit);
    return is;
  }
  grid = gd;
  return is;
}

std::ostream &operator<<(std::ostream &os, const Sudoku::Grid::Data &gd) {
  for (int i = 0; i < 81; ++i) {
    char c = '0' + gd.cells[i];
    if (!os.put(c))
      return os;
  }
  return os;
}

std::istream &operator>>(std::istream &is, Sudoku::Grid::Data &gd) {
  char buf[81];
  char c;

  do if (!is.get(c)) return is;
  while (isspace(c));
  is.unget();

  for (int i = 0; i < 81; ++i) {
    if (!is.get(c))
      return is;
    if (!std::isdigit(c)) {
      is.unget();
      is.setstate(std::ios::failbit);
      return is;
    }
    buf[i] = c;
  }

  for (int i = 0; i < 81; ++i) {
    gd.cells[i] = buf[i] - '0';
  }

  return is;
}

}; // namespace Sudoku


/**
 * @brief process sudokus passed to stdin
 *
 * examplary:
 * wget -q -O - http://school.maths.uwa.edu.au/~gordon/sudoku17 | ./sudoku --solve
 *
 * according to given example above, if piped from local for max. throughput
 * (and output piped to /dev/null), it takes less than 10 seconds to scan and solve
 * approx. 50.000 sudokus with minimum number of 17 fields given.
 * derived on a 1.8Ghz machine.
 * note that for the example above all given sudokus have been solved twice,
 * one times for scanning (implies validation, thus solving)
 * and one times to derive its one permutation (uniqueness expected).
 * so the solver logic itself is able to solve one sudoku within 100 microseconds.
 *
 * creating sudokus examplary:
 * ./sudoku --random 1 | ./sudoku --randet 10000 | ./sudoku --maxcl 22 |
 *      ./sudoku --friends -1 | ./sudoku --sort | ./sudoku --friends -1 | ./sudoku --sort |
 *      ./sudoku --friends -1 | ./sudoku --sort | ./sudoku --friends -1 | ./sudoku --sort |
 *      ./sudoku --stat
 *
 * according to given sample above one random sudoku gets generated,
 * then 10000 random determinations get derived from (produces sudokus with clues in range 21-29 normally),
 * then all sudokus with more then 22 clues get cut
 * and finally friend variations with one clue less get derived 4 times in sequence.
 * the result will be sudokus with clues in range 19-20 normally.
 */
int main(int argc, const char * const *argv) {
  if ((argc > 1) && (std::strcmp(argv[1], "--help") == 0)) {
    std::cout << *argv << " {--solve | --norm | --print | --random [N] | --randet [N] | --friends [-1 | only] | --sort | --stat | --maxcl N}" << std::endl;
  } else if ((argc > 1) && (std::strcmp(argv[1], "--solve") == 0)) {    // solve sudokus
    Sudoku::Grid grid;
    while (std::cin >> grid) {
      __int64_t N;
      grid.Permutate(N,
          [](__int64_t N, const Sudoku::Grid &grid) -> bool {std::cout << grid << std::endl; return N == 1;});
      if (N != 1) {
        std::cerr << grid << " not unique!" << std::endl;
      }
    }
  } else if ((argc > 1) && (std::strcmp(argv[1], "--norm") == 0)) {     // normalize sudokus
    Sudoku::Grid::Data gd, gn;
    while (std::cin >> gd) {
      std::cout << gd.Normalize(gn) << std::endl;
    }
  } else if ((argc > 1) && (std::strcmp(argv[1], "--print") == 0)) {    // print sudokus formatted
    Sudoku::Grid::Data gd;
    bool next = false;
    while (std::cin >> gd) {
      if (next)
        std::cout << std::endl << std::endl;
      gd.Print(std::cout);
      next = true;
    }
  } else if ((argc > 1) && (std::strcmp(argv[1], "--random") == 0)) {   // generate random sudokus
    Sudoku::Grid grid;
    __int64_t N, Nmax = 1000;
    if (argc > 2) {
      Nmax = 0;
      const char *q = argv[2];
      int c;
      while ((c = *q++) != '\0') {
        if (!std::isdigit(c)) {
          std::cerr << "option failure." << std::endl;
          return 1;
        }
        Nmax *= 10;
        Nmax += (c - '0');
      }
    }
    grid.Random();
    grid.Permutate(N,
        [&](__int64_t n, const Sudoku::Grid &g) -> bool {std::cout << g << std::endl; return n < Nmax;});
  } else if ((argc > 1) && (std::strcmp(argv[1], "--friends") == 0)) {  // find sudoku "friends"
    Sudoku::Grid grid, Gmax;
    std::set<Sudoku::Grid::Data> Gmin;
    bool less = ((argc > 2) && (std::strcmp(argv[2], "-1") == 0));      // only real subset friends (with less clues)
    bool only = ((argc > 2) && (std::strcmp(argv[2], "only") == 0));    // only friends, not origin
    int cmin = 81, cless = 81;
    while (std::cin >> grid) {
      __int64_t N;
      grid.Permutate(N,
          [&](__int64_t, const Sudoku::Grid &g) -> bool {Gmax = g; return false;});
      int clues = 0;
      bool haveless = false;
      for (int i = 0; i < 81; ++i) {
        Sudoku::Cell *cell = grid[i];
        if (!cell->IsFixed())
          continue;
        ++clues;
        unsigned bit = cell->Unfix();
        if (less && grid.IsUnique()) {
          haveless = true;
          std::cout << grid << std::endl;
        } else for (int j = 0; j < 81; ++j) {
          if (j == i)
            continue;
          Sudoku::Cell *cell = grid[j];
          if (cell->IsFixed())
            continue;
          cell->Fix(Gmax[j]->GetFixbit());
          if (grid.IsUnique()) {
            if (less) for (int k = 0; k < 81; ++k) {
              if (k == j)
                continue;
              Sudoku::Cell *cell = grid[k];
              if (!cell->IsFixed())
                continue;
              unsigned bit = cell->Unfix();
              if (grid.IsUnique()) {
                haveless = true;
                std::cout << grid << std::endl;
              }
              cell->Fix(bit);
            } else {
              std::cout << grid << std::endl;
            }
          }
          cell->Unfix();
        }
        cell->Fix(bit);
      }
      if (less) {
        if (clues < cmin) {
          cmin = clues;
          Gmin.clear();
        }
        if ((clues - (haveless ? 1 : 0)) < cless) {
          cless = clues - (haveless ? 1 : 0);
          if (cless < cmin)
            Gmin.clear();
        }
        if ((clues == cmin) && !(cless < cmin))
          Gmin.insert(Sudoku::Grid::Data(grid));
      } else if (!only) {
        std::cout << grid << std::endl;
      }
    }
    if (less && !(cless < cmin)) {
      // dump origin set, for which we have not found friends with less clues.
      for (std::set<Sudoku::Grid::Data>::const_iterator it = Gmin.begin(); it != Gmin.end(); ++it) {
        std::cout << *it << std::endl;
      }
    }
  } else if ((argc > 1) && (std::strcmp(argv[1], "--sort") == 0)) {     // sort sudokus
    Sudoku::Grid::Data gd;
    std::set<Sudoku::Grid::Data> gset;
    while (std::cin >> gd) {
      if (gset.find(gd) == gset.end()) {
        std::cout << gd << std::endl;
        gset.insert(gd);
      }
    }
  } else if ((argc > 1) && (std::strcmp(argv[1], "--stat") == 0)) {     // generate stats to stderr
    Sudoku::Grid::Data gd;
    __int64_t N = 0;
    int cmax = 0, cmin = 81;
    while (std::cin >> gd) {
      ++N;
      int i, clues;
      for (i = 0, clues = 0; i < 81; ++i) {
        if (gd[i] != 0)
          ++clues;
      }
      if (clues > cmax)
        cmax = clues;
      if (clues < cmin)
        cmin = clues;
      std::cout << gd << std::endl;
    }
    std::cerr << "# count = " << N;
    if (N != 0) {
      std::cerr << ", clues = " << cmin;
      if (cmax != cmin)
        std::cerr << "-" << cmax;
    }
    std::cerr << std::endl;
  } else if ((argc > 1) && (std::strcmp(argv[1], "--randet") == 0)) {   // random determination
    Sudoku::Grid grid, g;
    int N = 1;
    if (argc > 2) {
      N = 0;
      const char *q = argv[2];
      int c;
      while (((c = *q) != '\0') && std::isdigit(c)) {
        ++q;
        N *= 10;
        N += (c - '0');
      }
      if ((*q != '\0') || (N == 0)) {
        std::cerr << "option failure." << std::endl;
        return 1;
      }
    }
    std::srand(std::time(nullptr));
    while (std::cin >> grid) {
      for (int n = N; --n >= 0;) {
        g = grid;
        int idx[81];
        for (int i = 0; i < 81; ++i)
          idx[i] = i;
        for (int k = 6561; --k >= 0;) {
          int i0 = std::rand() % 81;
          int i1 = std::rand() % 81;
          if (i0 != i1) {
            int i = idx[i0];
            idx[i0] = idx[i1];
            idx[i1] = i;
          } else {
            ++k;
          }
        }
        for (int i = 0; i < 81; ++i) {
          Sudoku::Cell *cell = g[idx[i]];
          unsigned bit = cell->Unfix();
          if ((bit == 0) || g.IsUnique())
            continue;
          cell->Fix(bit);
        }
        std::cout << g << std::endl;
      }
    }
  } else if ((argc > 1) && (std::strcmp(argv[1], "--maxcl") == 0)) {    // filter sudokus (with more than max. clues given)
    Sudoku::Grid::Data gd;
    int N = 0;
    const char *q = (argc > 2) ? argv[2] : nullptr;
    if (q != nullptr) {
      int c;
      while (((c = *q) != '\0') && std::isdigit(c)) {
        ++q;
        N *= 10;
        N += (c - '0');
      }
    }
    if ((N == 0) || (*q != '\0')) {
      std::cerr << "option failure." << std::endl;
      return 1;
    }
    while (std::cin >> gd) {
      int i, clues;
      for (i = 0, clues = 0; i < 81; ++i) {
        if (gd[i] != 0)
          ++clues;
      }
      if (clues <= N)
        std::cout << gd << std::endl;
    }
  } else {
    std::cerr << "option failure." << std::endl;
    return 1;
  }

  return 0;
}
