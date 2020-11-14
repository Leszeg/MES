# Opis do lab :
# Całkując wykorzystuje zmienne w układzie lokalnym (ksi i eta)
# Funkcje dla jakich liczę całki wpisywane w funkcji integral
# Obiekt klasy element ma dodaną tablicę nodes w której na sztywno podałem współrzędne ksi i eta
# Nie ma powiązania z wcześniejszymi labolatoriami jeżeli chodzi o węzły.


    # Opis:
    # FEM_Grid - zawiera listę elementów oraz węzłów
    #     Element zna współrzędne swoich punktów i ich ID
    #     Węzły nie wiedzą do jakiego elementu należą i nie znają swojego ID
    #     Klasa FEM_Grid rysuje siatkę na podstawie swojej listy elementów.
    #     Lista w klasie Element jest używana tylko do obliczeń
    # GlobalData - dane odczytane z pliku
    # Element - Zawiera
    #     tablice z ID swoich węzłów(nodes_ID)
    #     tablice z współrzędnymi swoich węzłów(nodes)
    #     zmienną informującą o ilości węzłów w elemencie
    # Node - zawiera współrzędne punktów
    # Plot - drukuje wykres FEM_Grid
    # Read - czyta dane z pliku

    # Lab 2 - całkowanie numeryczne
    # Całkuje element 4 węzłowy