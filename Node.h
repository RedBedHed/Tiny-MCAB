//
// Created by evcmo on 2/6/2022.
//

#ifndef UNTITLED3_NODE_H
#define UNTITLED3_NODE_H

#include <vector>

namespace bit {

    /**
     * The Alliances, enumerated.
     */
    enum Alliance : uint8_t
    { X, O };

    struct Node final {
        int move = -1;
        double v = 0;
        double n = 0;
        double d = 0;
        Alliance a = X;
        Node *parent = nullptr;
        std::vector<Node *> x; 
        // Not space efficient, 
        // but who really cares?
        inline Node(
        const int m, 
        const Alliance a, 
        Node* const p) : 
        move(m), 
        a(a), 
        parent(p) { }
    };
}

#endif //UNTITLED3_NODE_H
