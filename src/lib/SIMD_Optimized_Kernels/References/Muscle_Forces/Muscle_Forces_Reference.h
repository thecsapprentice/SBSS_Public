template<class T>
void Muscle_Forces_Reference(T f[3][8], const T fiber[3],
                             const T Ffiber[3], const T c1,
                             const T one_over_h,
                             const T cell_volume);

template<class T>
bool Muscle_Forces_Compare(const T f[3][8], const T f_reference[3][8]);
