import sys
import wxgen.driver
import wxgen.driver_db


def main():
   wxgen.driver.run(sys.argv)

def main_db():
   wxgen.driver_db.run(sys.argv)

if __name__ == '__main__':
   main()
