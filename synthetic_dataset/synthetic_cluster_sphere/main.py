#!/usr/bin/python3
import sys
import numpy as np
import matplotlib.pyplot as plt


def save_list(file_name, L):
    f = open(file_name + ".txt", "w")
    f.write(",".join([str(el) for el in L]))
    f.close()


def save_matrix(file_name, M):
    f = open(file_name + ".txt", "w")
    f.write(
        "\n".join(
            [",".join([str(M[row][col]) for col in M[row].keys()]) for row in M.keys()]
        )
    )
    f.close()


def save_data(F, C, dAtoC, dFtoF):
    save_list("F", F)
    save_list("C", C)
    save_matrix("dAtoC", dAtoC)
    save_matrix("dFtoF", dFtoF)


def save_text(file_name, text):
    f = open(file_name + ".txt", "w")
    f.write(text)
    f.close()


def save_average_distances(A, B, C, dAtoC, dFtoF):
    text = "AA, BB, AB, AC, BC\naverage:\n"
    average_distances = [0, 0, 0, 0, 0]
    min_distances = [
        dFtoF[A[0]][A[1]],
        dFtoF[B[0]][B[1]],
        dFtoF[B[0]][A[0]],
        dAtoC[A[0]][C[0]],
        dAtoC[B[0]][C[0]],
    ]

    max_distances = [
        dFtoF[A[0]][A[1]],
        dFtoF[B[0]][B[1]],
        dFtoF[B[0]][A[0]],
        dAtoC[A[0]][C[0]],
        dAtoC[B[0]][C[0]],
    ]
    for a1 in A:
        for a2 in A:
            average_distances[0] += dFtoF[a1][a2]
            if a1 != a2:
                min_distances[0] = min(min_distances[0], dFtoF[a1][a2])
                max_distances[0] = max(max_distances[0], dFtoF[a1][a2])
    average_distances[0] = average_distances[0] / (2 * (len(A) - 1) * len(A))
    for b1 in B:
        for b2 in B:
            average_distances[1] += dFtoF[b1][b2]
            if b1 != b2:
                min_distances[1] = min(min_distances[1], dFtoF[b1][b2])
                max_distances[1] = max(max_distances[1], dFtoF[b1][b2])
    average_distances[1] = average_distances[1] / (2 * (len(B) - 1) * len(B))
    for b in B:
        for a in A:
            average_distances[2] += dFtoF[a][b]
            min_distances[2] = min(min_distances[2], dFtoF[a][b])
            max_distances[2] = max(max_distances[2], dFtoF[a][b])
    average_distances[2] = average_distances[2] / (len(B) * len(A))
    for a in A:
        for c in C:
            average_distances[3] += dAtoC[a][c]
            min_distances[3] = min(min_distances[3], dAtoC[a][c])
            max_distances[3] = max(max_distances[3], dAtoC[a][c])
    average_distances[3] = average_distances[3] / (len(A) * len(C))
    for b in B:
        for c in C:
            average_distances[4] += dAtoC[b][c]
            min_distances[4] = min(min_distances[4], dAtoC[b][c])
            max_distances[4] = max(max_distances[4], dAtoC[b][c])
    average_distances[4] = average_distances[4] / (len(B) * len(C))
    text += ",".join([str(el) for el in average_distances])
    text += "\nmin:\n"
    text += ",".join([str(el) for el in min_distances])
    text += "\nmax:\n"
    text += ",".join([str(el) for el in max_distances])
    save_text("average_distances", text)


def create_on_sphere(offset, r, num, d):
    loc = []
    vectors = np.random.normal(0, 1, num * d)
    for i in range(num):
        # loc.append([r*np.cos(dir1[i]) + offsetX, r*np.sin(dir1[i]) + offsetY])
        coord = np.array(vectors[i * d : (i + 1) * d])
        coord = coord / np.linalg.norm(coord)
        loc.append(r * coord + offset)
    return loc


def create_in_sphere(offset, r, num, d):
    loc = []
    vectors = np.random.normal(0, 1, num * d)
    scalar = np.random.uniform(0, 1, num)
    for i in range(num):
        # loc.append([r*np.cos(dir1[i]) + offsetX, r*np.sin(dir1[i]) + offsetY])
        coord = np.array(vectors[i * d : (i + 1) * d])
        coord = coord / np.linalg.norm(coord)
        loc.append(r * coord * scalar[i] + offset)
    return loc


def euclidean_norm(p1, p2):
    p1 = np.array(p1)
    p2 = np.array(p2)
    return np.linalg.norm(p1 - p2)


def create_dAtoC(F, C):
    dAtoC = {a: {} for a in range(len(F) + len(C))}
    for (f, f_loc) in F:
        for (c, c_loc) in C:
            dAtoC[f][c] = euclidean_norm(f_loc, c_loc)
    for (c1, c1_loc) in C:
        for (c2, c2_loc) in C:
            dAtoC[c1][c2] = euclidean_norm(c1_loc, c2_loc)
    return dAtoC


def create_dFtoF(F):
    dFtoF = {f: {} for f in range(len(F))}
    for (f1, f1_loc) in F:
        for (f2, f2_loc) in F:
            dFtoF[f1][f2] = euclidean_norm(f1_loc, f2_loc)
    return dFtoF


def create_bad_loc(bad_r, num_bad, d):
    offset = np.zeros(d)
    bad_offsets = create_on_sphere(offset, bad_r, len(num_bad), d)
    bad_loc = []
    for i in range(len(num_bad)):
        bad_loc = bad_loc + create_in_sphere(np.array(bad_offsets[i]), 1, num_bad[i], d)
    return bad_loc


def create_synthetic_cluster_sphere(
    opt_r, client_r, bad_r, num_good, num_clients, num_bad, d
):
    good_offset = np.zeros(d)
    client_offset = np.zeros(d)
    good_loc = create_in_sphere(good_offset, opt_r, num_good, d)
    client_loc = create_in_sphere(client_offset, client_r, num_clients, d)
    bad_loc = create_bad_loc(bad_r, num_bad, d)
    print("Created locations")

    F = [[i, good_loc[i]] for i in range(len(good_loc))]
    F = F + [[i + len(good_loc), bad_loc[i]] for i in range(len(bad_loc))]
    C = [
        [i + len(good_loc) + len(bad_loc), client_loc[i]]
        for i in range(len(client_loc))
    ]
    print("Got F and C")

    dAtoC = create_dAtoC(F, C)
    dFtoF = create_dFtoF(F)
    print("Calculated distances")
    plt.gca().set_aspect("equal", adjustable="box")
    plt.scatter([el[0] for el in bad_loc], [el[1] for el in bad_loc], color="red")
    plt.scatter(
        [el[0] for el in client_loc], [el[1] for el in client_loc], color="black"
    )
    plt.scatter([el[0] for el in good_loc], [el[1] for el in good_loc], color="green")
    plt.axis("square")
    plt.savefig("sphere_plot.png")
    F = [el[0] for el in F]
    C = [el[0] for el in C]
    save_data(F, C, dAtoC, dFtoF)

    A = [i for i in range(num_good)]
    B = [i for i in range(num_good, num_bad[0])]
    save_average_distances(A, B, C, dAtoC, dFtoF)


def main(arguments):
    (opt_r, client_r, bad_r, num_good, num_clients, num_bad, d) = arguments
    create_synthetic_cluster_sphere(
        opt_r, client_r, bad_r, num_good, num_clients, num_bad, d
    )


def ready_argv(argv):
    argv = sys.argv[1:]
    if len(argv) != 7:
        print("Not enough arguments")
    input_num = []
    for i in range(0, 3):
        input_num.append(float(argv[i]))
    for i in range(3, 5):
        input_num.append(int(argv[i]))
    input_num.append([int(el) for el in argv[5].split(",")])
    input_num.append(int(argv[6]))
    return input_num


if __name__ == "__main__":
    main(ready_argv(sys.argv))
