#include <iostream>

bool is_palindrome(std::string text) {
  
    std::string reverse = "";
  
    for (int i = text.size() - 1; i >= 0; i--) {
        reverse += text[i];
    } 

    if (reverse == text) {
        return true;
    } else {
        return false;
    }

}

int main() {
  
  std::cout << is_palindrome("madam") << "\n";
  std::cout << is_palindrome("ada") << "\n";
  std::cout << is_palindrome("lovelace") << "\n";
  
}