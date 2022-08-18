/// Cliff Sparks 2021
/// Game of Hex
/// C++ for C Programmers, Part B: Week 4 Homework 2

#include <iostream>
#include <string>
#include <random>
#include <exception>
#include <cassert>
#include <ctime>
//#include <cstdlib>

using namespace std;

/// === hsize is the board size === ///
const int hsize = 11;   // (now set at compile-time for ~10% optimization)
/// === maxN is the number of Monte-Carlo iterations === ///
const int maxN = 1000;   // (higher for better but slower play)
// (maxN = 1000 takes about 17 seconds per move on a Core i5 11th gen laptop
// with no (GNU) compiler optimizations, and about 5.5 seconds per move with
// several (GNU) compiler optimizations enabled.

const int maxhsize = 11;   // max board size
const int hnodes = hsize*hsize + 2;   // graph size

enum class Colour:short {None = -1, Blue, Red, Both};
Colour Blue = Colour::Blue;
Colour Red = Colour::Red;
short b = static_cast<short>(Blue);
short r = static_cast<short>(Red);

// C++ random
default_random_engine generator(time(nullptr));
uniform_real_distribution<double> distribution(0.0, 1.0);

class HexBoard {
  public:
    HexBoard();   // default constructor
    HexBoard(const HexBoard&);   // copy constructor
    ~HexBoard();   // destructor

    int getSize() const;   // size of board
    int getNodes() const;   // number of nodes
    // check whether a piece is on the board
    bool isPiece(const Colour, const int i, const int j) const;
    // returns a connection from the graph
    bool isCon(const Colour, const int m, const int n) const;
    // piece to output to the screen
    string screenPiece (const bool, const bool) const;
    // overload << operator to display board on the screen
    friend ostream& operator<< (ostream&, const HexBoard&);
    void printCon(const Colour) const;   // display the connectivity graph
    // make connections in the graph based on the new piece position (i, j)
    void makeCon(const Colour, const int i, const int j);
    void userPiece(const Colour);   // user puts a piece on the board
    // Depth First Search recursive algorithm
    void depthFS(const Colour, bool*, int*, const int, const int) const;
    bool isWon(const Colour) const;   // has a colour won
    bool displayWin() const;   // checks for and displays the winning colour

    // AI playing section...
    inline double prob() const;   // random number in the interval [0, 1)
    int remSquares(int*, const int) const;   // empty squares remaining
    // computer makes a random move in the remaining empty 'squares'
    bool randPiece(int index[], const int, const Colour);
    // copy pieces from another board (connection graph not copied)
    void copyPieces(const HexBoard& board);
    void reCon();   // re-calculate connections
    void MonteCarlo(const Colour);   // Monte Carlo simulations

  private:
    bool*** hpiece;  // blue/red piece arrays
    bool*** hcon;  // blue/red connectivity large arrays
};

// Default constructor
HexBoard::HexBoard() {
    // Assertion will fail if board size is out of range
    assert(hsize >= 1 && hsize <= maxhsize);

    hpiece = new bool**[2];
    hpiece[b] = new bool*[hsize];   // blue
    hpiece[r] = new bool*[hsize];   // red
    hcon = new bool**[2];
    hcon[b] = new bool*[hnodes];   // blue
    hcon[r] = new bool*[hnodes];   // red

    for (int i = 0; i < hsize; ++i) {
        hpiece[b][i] = new bool[hsize];
        hpiece[r][i] = new bool[hsize];
        // Initialize default values
        for (int j = 0; j < hsize; ++j) {
            hpiece[b][i][j] = false;
            hpiece[r][i][j] = false;
        }
    }

    // using (m, n) for the connectivity graph
    // to distinguish it from the board (i, j)
    for (int m = 0; m < hnodes; ++m) {
        hcon[b][m] = new bool[hnodes];
        hcon[r][m] = new bool[hnodes];
        // Initialize default values
        for (int n = 0; n < hnodes; ++n) {
            hcon[b][m][n] = false;
            hcon[r][m][n] = false;
        }
    }

    //cout << "HexBoard default constructor called" << endl;
};

// Construct a HexBoard with data from another HexBoard
HexBoard::HexBoard(const HexBoard& board) {
    hpiece = new bool**[2];
    hpiece[b] = new bool*[hsize];   // blue
    hpiece[r] = new bool*[hsize];   // red
    hcon = new bool**[2];
    hcon[b] = new bool*[hnodes];   // blue
    hcon[r] = new bool*[hnodes];   // red

    for (int i = 0; i < hsize; ++i) {
        hpiece[b][i] = new bool[hsize];
        hpiece[r][i] = new bool[hsize];
        // Initialize default values
        for (int j = 0; j < hsize; ++j) {
            hpiece[b][i][j] = board.isPiece(Blue, i, j);
            hpiece[r][i][j] = board.isPiece(Red, i, j);
        }
    }

    // using (m, n) for the connectivity graph
    // to distinguish it from the board (i, j)
    for (int m = 0; m < hnodes; ++m) {
        hcon[b][m] = new bool[hnodes];
        hcon[r][m] = new bool[hnodes];
        // Initialize default values
        for (int n = 0; n < hnodes; ++n) {
            hcon[b][m][n] = board.isCon(Blue, m, n);
            hcon[r][m][n] = board.isCon(Red, m, n);
        }
    }

    //cout << "HexBoard copy constructor called" << endl;
};

// Destructor
HexBoard::~HexBoard() {
    for (int i = 0; i < hsize; i++) {
        delete[] hpiece[b][i];
        delete[] hpiece[r][i];
    }
    delete[] hpiece[b];
    delete[] hpiece[r];
    delete[] hpiece;

    for (int m = 0; m < hnodes; m++) {
        delete[] hcon[b][m];
        delete[] hcon[r][m];
    }
    delete[] hcon[b];
    delete[] hcon[r];
    delete[] hcon;

    //cout << "HexBoard destructor called" << endl;
};

// Get the size of the board
int HexBoard::getSize() const {
    return hsize;
}

// Get the nodes in the graph
int HexBoard::getNodes() const {
    return hnodes;
}

// Check whether a piece is on the board
bool HexBoard::isPiece(const Colour Col, const int i, const int j) const {
    short c = static_cast<short>(Col);

    if ((Col != Blue && Col != Red)
            || i < 0 || i >= hsize || j < 0 || j >= hsize) {
        cout << "Piece coordinates out of range!" << endl;
        return false;
    }
    else
        return hpiece[c][i][j];
}

// Returns true if connection exists in the graph
bool HexBoard::isCon(const Colour Col, const int m, const int n) const {
    short c = static_cast<short>(Col);

    if ((Col != Blue && Col != Red)
            || m < 0 || m >= hnodes || n < 0 || n >= hnodes) {
        cout << "Connection coordinates out of range!" << endl;
        return false;
    }
    else
        return hcon[c][m][n];
}

// Piece to output to the screen
string HexBoard::screenPiece(const bool bl, const bool rd) const {
    string piece;

    if (!bl && !rd)
        piece = ".";   // no piece
    else if (bl && !rd)
        piece = "B";   // blue piece
    else if (!bl && rd)
        piece = "R";   // red piece
    else   // (bl && rd)
        piece = "#";   // this should not happen

    return piece;
}

// Overload << operator to display board on the screen
ostream& operator<< (ostream& out, const HexBoard& board) {
    const int s = hsize;
    string offset = "";

    for (int i = 0; i < s - 1; ++i) {
        out << offset;
        for (int j = 0; j < s - 1; ++j) {
            out << board.screenPiece(board.hpiece[b][i][j],
                                        board.hpiece[r][i][j]);
            out << " - ";
        }
        out << board.screenPiece(board.hpiece[b][i][s - 1],
                                    board.hpiece[r][i][s - 1]);
        out << endl << offset << " ";

        for (int j = 0; j < s - 1; ++j) {
            out << "\\ / ";   // double backslash to override escape sequence
        }
        out << "\\" << endl;   // double backslash to override escape sequence
        offset += "  ";
    }

    out << offset;
    for (int j = 0; j < s - 1; ++j) {
        out << board.screenPiece(board.hpiece[b][s - 1][j],
                                    board.hpiece[r][s - 1][j]);
        out << " - ";
    }
    out << board.screenPiece(board.hpiece[b][s - 1][s - 1],
                                board.hpiece[r][s - 1][s - 1]);

    return out;
}

// Print the connectivity graph
void HexBoard::printCon(const Colour Col) const {
    short c = static_cast<short>(Col);
    for (int m = 0; m < hnodes; m++) {
        for (int n = 0; n < hnodes; n++)
            cout << hcon[c][m][n] << " ";
        cout << endl;
    }
    cout << endl;
}

// Make connections in the graph based on the new piece position (i, j)
void HexBoard::makeCon(const Colour Col, const int i, const int j) {
    short c = static_cast<short>(Col);

    if (Col == Red) {
        // Graph node 0 represents the top of the board
        // Let 'North' be the direction from the bottom to the top of the board
        if (i == 0) {
            hcon[c][0][j + 1] = true;
            hcon[c][j + 1][0] = true;   // undirected graph...
        }

        // Look for a 'North-West' piece
        if (i - 1 != -1 && hpiece[c][i - 1][j]) {
            hcon[c][(i - 1)*hsize + j + 1][i*hsize + j + 1] = true;
            hcon[c][i*hsize + j + 1][(i - 1)*hsize + j + 1] = true;
        }
        // Look for a 'North-East' piece
        if ((i - 1 != -1 && j + 1 != hsize) && hpiece[c][i - 1][j + 1]) {
            hcon[c][(i - 1)*hsize + j + 2][i*hsize + j + 1] = true;
            hcon[c][i*hsize + j + 1][(i - 1)*hsize + j + 2] = true;
        }
        // Look for a 'West' piece
        if (j - 1 != -1 && hpiece[c][i][j - 1]) {
            hcon[c][i*hsize + j][i*hsize + j + 1] = true;
            hcon[c][i*hsize + j + 1][i*hsize + j] = true;
        }
        // Look for an 'East' piece
        if (j + 1 != hsize && hpiece[c][i][j + 1]) {
            hcon[c][i*hsize + j + 2][i*hsize + j + 1] = true;
            hcon[c][i*hsize + j + 1][i*hsize + j + 2] = true;
        }
        // Look for a 'South-West' piece
        if ((i + 1 != hsize && j - 1 != -1) && hpiece[c][i + 1][j - 1]) {
            hcon[c][(i + 1)*hsize + j][i*hsize + j + 1] = true;
            hcon[c][i*hsize + j + 1][(i + 1)*hsize + j] = true;
        }
        // Look for a 'South-East' piece
        if (i + 1 != hsize && hpiece[c][i + 1][j]) {
            hcon[c][(i + 1)*hsize + j + 1][i*hsize + j + 1] = true;
            hcon[c][i*hsize + j + 1][(i + 1)*hsize + j + 1] = true;
        }

        // Graph node (hnodes - 1) represents the bottom of the board
        if (i == hsize - 1) {
            hcon[c][hnodes - 1][i*hsize + j + 1] = true;
            hcon[c][i*hsize + j + 1][hnodes - 1] = true;
        }
    }
    else {   // Col == Blue
        // Graph node 0 represents the left of the board
        // Let 'North' be the direction from the left to the right of the board
        if (j == 0) {
            hcon[c][0][i + 1] = true;
            hcon[c][i + 1][0] = true;   // undirected graph...
        }

        // Look for a 'South' piece
        if (j - 1 != -1 && hpiece[c][i][j - 1]) {
            hcon[c][(j - 1)*hsize + i + 1][j*hsize + i + 1] = true;
            hcon[c][j*hsize + i + 1][(j - 1)*hsize + i + 1] = true;
        }
        // Look for a 'South-East' piece
        if ((j - 1 != -1 && i + 1 != hsize) && hpiece[c][i + 1][j - 1]) {
            hcon[c][(j - 1)*hsize + i + 2][j*hsize + i + 1] = true;
            hcon[c][j*hsize + i + 1][(j - 1)*hsize + i + 2] = true;
        }
        // Look for a 'South-West' piece
        if (i - 1 != -1 && hpiece[c][i - 1][j]) {
            hcon[c][j*hsize + i][j*hsize + i + 1] = true;
            hcon[c][j*hsize + i + 1][j*hsize + i] = true;
        }
        // Look for a 'North-East' piece
        if (i + 1 != hsize && hpiece[c][i + 1][j]) {
            hcon[c][j*hsize + i + 2][j*hsize + i + 1] = true;
            hcon[c][j*hsize + i + 1][j*hsize + i + 2] = true;
        }
        // Look for a 'North-West' piece
        if ((j + 1 != hsize && i - 1 != -1) && hpiece[c][i - 1][j + 1]) {
            hcon[c][(j + 1)*hsize + i][j*hsize + i + 1] = true;
            hcon[c][j*hsize + i + 1][(j + 1)*hsize + i] = true;
        }
        // Look for a 'North' piece
        if (j + 1 != hsize && hpiece[c][i][j + 1]) {
            hcon[c][(j + 1)*hsize + i + 1][j*hsize + i + 1] = true;
            hcon[c][j*hsize + i + 1][(j + 1)*hsize + i + 1] = true;
        }

        // Graph node (hnodes - 1) represents the right of the board
        if (j == hsize - 1) {
            hcon[c][hnodes - 1][j*hsize + i + 1] = true;
            hcon[c][j*hsize + i + 1][hnodes - 1] = true;
        }
    }
}

// User puts a piece on the board
void HexBoard::userPiece(const Colour Col) {
    short c = static_cast<short>(Col);
    string scol[2] = {"Blue", "Red"};
    string si, sj;

    while (true) {
        cout << "Enter row for your " << scol[c] << " piece: ";
        cin >> si;
        cout << "Enter column for your " << scol[c] << " piece: ";
        cin >> sj;

        int i = -1, j = -1;
        try {
            i = stoi(si);
            j = stoi(sj);
        }
        catch(...) {};

        if (i < 1 || i > hsize || j < 1 || j > hsize)
            cout << "Move out of range!" << endl;
        else if (hpiece[1 - c][i - 1][j - 1])
            cout << "There is a " << scol[1 - c] << " piece there!" << endl;
        else if (hpiece[c][i - 1][j - 1])
            cout << "There is already a " << scol[c] << " piece there!" << endl;
        else {
            hpiece[c][i - 1][j - 1] = true;
            makeCon(Col, i - 1, j - 1);
            cout << "Your " << scol[c] << " piece has been placed at ("
                << i << ", " << j << ")." << endl;
            break;
        }
    }
}

// Depth First Search (DFS) recursive algorithm
void HexBoard::depthFS(const Colour Col, bool connected[], int index[],
                        const int m, const int id) const {
    short c = static_cast<short>(Col);
    connected[m] = true;
    index[m] = id;
    for (int n = 0; n < hnodes; n++) {
        if (!hcon[c][m][n]) continue;
        else if (!connected[n])
            depthFS(Col, connected, index, n, id);
    }
}

// Has a colour won (all board sizes)
bool HexBoard::isWon(const Colour Col) const {
    short c = static_cast<short>(Col);
    bool connected[hnodes];
    int index[hnodes];
    for (int k = 0; k < hnodes; ++k) {
        connected[k] = false;
        index[k] = -1;
    }

    int count = 0;
    for (int k = 0; k < hnodes; k++) {
        if (!connected[k]) {
            depthFS(Col, connected, index, k, count + 1);
            count++;
            // cout << "Subgraph count: " << count << endl;
        }
    }

    // true if the top(/left) of the board is connected to the bottom(/right)
    return (index[0] != -1 && index[0] == index[hnodes - 1]);
}

// Checks for and displays the winning colour
bool HexBoard::displayWin() const {
    int won = false;
    if (isWon(Red)) {
        // ASCII code 7 is 'bel' - thus plays a sound when the game is won
        cout << "Red has won!" << static_cast<char>(7) << endl;
        won = true;
    }
    else if (isWon(Blue)) {
        cout << "Blue has won!" << static_cast<char>(7) << endl;
        won = true;
    }

    return won;
}

// Random number in the interval [0, 1)
inline double HexBoard::prob() const {
    // interval [0, 1]
    //return static_cast<double>(rand())/static_cast<double>(RAND_MAX);
    // interval [0, 1)
    return distribution(generator);
}

// Empty 'squares' remaining on the board
int HexBoard::remSquares(int index[], const int sqrtsize) const {
    int i, j, k;
    int hsquares = -1;   // hex squares remaining
    for (k = 0; k < sqrtsize*sqrtsize; ++k) {
        i = k / sqrtsize;
        j = k % sqrtsize;
        if (!hpiece[0][i][j] && !hpiece[1][i][j]) {
            ++hsquares;
            index[hsquares] = k;
        }
    }

    return hsquares;
}

// Computer makes a random move in the remaining empty 'squares'
bool HexBoard::randPiece(int index[], const int sqrtsize, const Colour Col) {
    short c = static_cast<short>(Col);
    int i, j, k;

    int hsquares = remSquares(index, sqrtsize);
    if (hsquares == -1) return false;

    k = static_cast<int>(prob() * static_cast<double>(hsquares + 1));
    // unlikely, only necessary for prob() interval [0, 1]
    //if (k == hsquares + 1) k = hsquares;
    i = index[k] / hsize;
    j = index[k] % hsize;
    hpiece[c][i][j] = true;

    return true;
}

// Copy pieces from another board (connection graph not copied)
void HexBoard::copyPieces(const HexBoard& board) {
    for (int i = 0; i < hsize; ++i) {
        for (int j = 0; j < hsize; ++j) {
            hpiece[b][i][j] = board.isPiece(Blue, i, j);
            hpiece[r][i][j] = board.isPiece(Red, i, j);
        }
    }
}

// Re-calculate connections
void HexBoard::reCon() {
    for (int m = 0; m < hnodes; ++m) {
        for (int n = 0; n < hnodes; ++n) {
            hcon[b][m][n] = false;
            hcon[r][m][n] = false;
        }
    }

    for (int i = 0; i < hsize; ++i) {
        for (int j = 0; j < hsize; ++j) {
            if (hpiece[b][i][j]) makeCon(Blue, i, j);
            if (hpiece[r][i][j]) makeCon(Red, i, j);
        }
    }
}

// Monte Carlo simulations
void HexBoard::MonteCarlo(const Colour Col) {
    short c = static_cast<short>(Col);
    int maxcol[2] = {0, 0};
    int index_one[hsize*hsize];
    int hsquares = remSquares(index_one, hsize);
    int besti = -1, bestj = -1;

    HexBoard startboard(*this);
    HexBoard monteboard(startboard);

    int index_two[hsize*hsize];
    for (int sq = 0; sq < hsquares + 1; ++sq) {
        int i = index_one[sq] / hsize;
        int j = index_one[sq] % hsize;
        startboard.hpiece[c][i][j] = true;
        // (set graph connections later)

        int ncol[2] = {0, 0};
        for (int N = 0; N < maxN; ++N) {
            int k = 1 - c;
            monteboard.copyPieces(startboard);
            for (int iter = 0; ; ++iter) {
                Colour Kol = static_cast<Colour>(k);
                if (!monteboard.randPiece(index_two, hsize, Kol)) break;
                k = 1 - k;
                //assert(iter != hsquares);
            }
            monteboard.reCon();   // now calculate connections
            if (monteboard.isWon(Red)) ++ncol[r];
            else ++ncol[b];   // if red hasn't won, blue must've won
        }

        if (ncol[c] > maxcol[c]) {
            maxcol[c] = ncol[c];
            besti = i;
            bestj = j;
        }

        startboard.hpiece[c][i][j] = false;   // remove piece
    }

    /*
    int redwon, bluewon;
    if (c == r) {
        redwon = maxcol[c];
        bluewon = maxN - redwon;
    }
    else {   // c == b
        bluewon = maxcol[c];
        redwon = maxN - bluewon;
    }
    cout << "Red won " << redwon << " times; Blue won "
            << bluewon << " times." << endl;
    */

    hpiece[c][besti][bestj] = true;
    makeCon(Col, besti, bestj);
    string scol[2] = {"Blue", "Red"};
    cout << "Computer put a " << scol[c] << " piece at (" << besti + 1
            << ", " << bestj + 1 << ")" << endl;
    cout << scol[c] << " has a " << (100 * maxcol[c])/maxN
            << "% chance of winning." << endl;
}

int main() {
    //srand(time(nullptr));   // seed rand() for HexBoard::prob() method

    cout << "======= Game of Hex =======" << endl;
    cout << "Blue path is left <-> right." << endl;
    cout << "Red path is top <-> bottom." << endl;
    cout << "Blue goes first." << endl;
    cout << "Board size is " << hsize << " by " << hsize << "." << endl;

    HexBoard board;
    cout << endl << board << endl << endl;

    char colour;
    cout << "What colour would you like to be? ('r' for Red, 'b' for Blue) ";
    cin >> colour;
    if (colour == 'r')
        cout << "You are Red!" << endl;
    else {
        colour = 'b';
        cout << "You are Blue!" << endl;
    }
    cout << endl;

    if (colour == 'r') {
        while (true) {
            cout << "Computer's move..." << endl;
            board.MonteCarlo(Blue);
            cout << endl << board << endl << endl;
            if (board.displayWin()) break;
            board.userPiece(Red);
            cout << endl << board << endl << endl;
            if (board.displayWin()) break;
        }
    }
    else {   // colour == 'b'
        while (true) {
            board.userPiece(Blue);
            cout << endl << board << endl << endl;
            if (board.displayWin()) break;
            cout << "Computer's move..." << endl;
            board.MonteCarlo(Red);
            cout << endl << board << endl << endl;
            if (board.displayWin()) break;
        }
    }

    return 0;
};
