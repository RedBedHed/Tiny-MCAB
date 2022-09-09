#include <iostream>
#include <ostream>
#include <ctime>
#include <iomanip>
#include "Board.h"
#include "Opponent.h"

using namespace bit;
using std::cout;
using std::cin;
using std::flush;

/**
 * 5-min Tic Tac Toe. Please don't judge ':D
 */
int main() {
    srand(time(NULL));
    Board b; char c; int i = -1; bool x, o;
    uint16_t hash = 0;
    do {
        cout << "\033[2J\033[H" << flush;
        cout << b;
        Node n(-1, O, nullptr);
        Node* r = &n;
        bool flag = true;
        while(!b.isFull() &&
              !(x = b.hasVictory<X>()) &&
              !(o = b.hasVictory<O>())) {
            cout << "\nMake a move (1-9)\n"
                 << "9 8 7\n"
                 << "6 5 4\n"
                 << "3 2 1\n>>_";
            while(!(cin >> i)) {
                cin.clear();
                cin.ignore(INT32_MAX,'\n');
            }
            if(i < 1 || i > 9 || b.occupiedSquare(i - 1))
                continue;    
            if(flag) {  
                b.mark<O>(i - 1);
                opponent::select<true>(&b, r);
            } else {
                r = opponent::child(&b, r, i - 1);
                b.mark<O>(i - 1);
            }
            if(b.isFull()) break;
            const long start = clock();
            const int m = opponent::search(&b, r, flag);
            b.mark<X>(m);
            const long end = clock() - start;
            cout << "\033[2J\033[H" << flush;
            cout << b;
            printf(
                    "I took %.6f seconds!",
                    ((double)end) / CLOCKS_PER_SEC
            );
            flag = false;
        }
        cout << '\n'
             << (x? "I win!":
                 o? "You win!":
                 "Tie!")
             << "\nplay again? (y/n)\n>>_";
        b.reset();
        hash = 0;
        opponent::destroyTree(&n);
        while(!(cin >> c)) {
            cin.clear();
            cin.ignore(INT32_MAX,'\n');
        }
    } while(c == 'y');
    return 0;
}
